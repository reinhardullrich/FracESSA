#include <fracessa/exact_candidate_solver.hpp>

#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace fracessa::search {

namespace {

size_t checked_reduced_entry_cache_size(size_t dimension)
{
    const size_t maximum = std::numeric_limits<size_t>::max();
    if (dimension != 0 && dimension > maximum / dimension)
        throw std::overflow_error("Reduced-entry cache size overflow");
    const size_t square = dimension * dimension;
    if (dimension != 0 && square > maximum / dimension)
        throw std::overflow_error("Reduced-entry cache size overflow");
    return square * dimension;
}

} // namespace

/*
 * Exact candidate solver for a symmetric payoff matrix.
 *
 * The parser has already cleared the rational denominators once for the whole game:
 *
 *     integer_game = game_denominator * game.
 *
 * The denominator is positive, so this scaling preserves every equality, inequality, and inertia sign used by candidate search and stability.
 * Each support can therefore stay in integer arithmetic until a successful candidate is written to the public rational result.
 */
exact_candidate_solver::exact_candidate_solver(const numeric::matrix_int& integer_game, const numeric::integer& game_denominator)
    : dimension_(integer_game.rows())
    , reduced_entry_cache_(checked_reduced_entry_cache_size(dimension_))
    , reduced_entry_cache_ready_(checked_reduced_entry_cache_size(dimension_), 0)
    , ffldlt_workspace_(dimension_)
    , integer_game_(integer_game)
    , game_denominator_(game_denominator)
{
    if (dimension_ > support::kMaxBitsetDimension) {
        support_indices_large_.reserve(dimension_);
        non_support_indices_large_.reserve(dimension_);
    }
}

void exact_candidate_solver::resize_reduced_system(size_t reduced_dimension)
{
    if (reduced_system_dimension_ == reduced_dimension) return;

    reduced_system_dimension_ = reduced_dimension;
    reduced_system_.resize(reduced_system_dimension_, reduced_system_dimension_);
    right_hand_side_.resize(reduced_system_dimension_, 1);
    solution_numerators_.resize(reduced_system_dimension_, 1);
}

/*
 * Eliminate the normalization and payoff border before factorization.
 *
 * Let m be the first strategy in support S and let Z have columns e_i-e_m for every other i in S. Every normalized support vector has the unique
 * form x=e_m+Z*y. Multiplying A_S*x=u*1 by Z^T eliminates u and gives
 *
 *     H*y = r,       H = Z^T*A_S*Z,       r = -Z^T*A_S*e_m.
 *
 * The stored integer game is d*A for one positive d. We therefore build d*H and d*r below. This integer system has the same solution and inertia
 * as the rational system, and H is nonsingular exactly when the original bordered candidate matrix is nonsingular.
 */
template<class Index>
void exact_candidate_solver::build_reduced_system(const Index* support_indices, size_t reduced_dimension)
{
    const size_t reference = support_indices[0];
    const auto reference_diagonal = integer_game_(reference, reference);
    const size_t reference_cache_offset = reference * dimension_ * dimension_;

    for (size_t row = 0; row < reduced_dimension; ++row) {
        const size_t i = support_indices[row + 1];
        auto row_offset = right_hand_side_(row, 0);
        row_offset.set_difference(reference_diagonal, integer_game_(i, reference));

        // The cache is physically dense so that fixing m and i once per row leaves only direct [j] access in the inner loop.
        // Zero is a valid reduced entry, hence the separate byte array records whether each exact value has been calculated.
        const size_t cache_row_offset = reference_cache_offset + i * dimension_;
        auto* cached_values = reduced_entry_cache_.data() + cache_row_offset;
        auto* cached_ready = reduced_entry_cache_ready_.data() + cache_row_offset;

        // The right-hand side d*r_i=A_mm-A_im is also the repeated row term in
        // d*H_ij=A_ij-A_mj+(A_mm-A_im). The symmetric fraction-free factorization reads only the lower triangle.
        for (size_t column = 0; column <= row; ++column) {
            const size_t j = support_indices[column + 1];
            auto& cached_value = cached_values[j];
            if (!cached_ready[j]) {
                cached_value.set_difference(integer_game_(i, j), integer_game_(reference, j));
                cached_value += row_offset;
                cached_ready[j] = 1;
            }
            reduced_system_(row, column) = cached_value;
        }
    }
}

/*
 * Return one entry of the integer-scaled reduced extended Hessian in difference coordinates:
 *
 *     d*(e_row-e_reference)^T*A*(e_column-e_reference)
 *       = dA_row,column-dA_reference,column+dA_reference,reference-dA_row,reference.
 *
 * The support Hessian H and the later G and Q blocks reuse these same entries, so calculate each triple only once.
 */
numeric::integer::const_reference exact_candidate_solver::reduced_entry(size_t reference, size_t row, size_t column)
{
    const size_t offset = reference * dimension_ * dimension_ + row * dimension_ + column;
    auto& value = reduced_entry_cache_[offset];
    if (!reduced_entry_cache_ready_[offset]) {
        value.set_difference(integer_game_(row, column), integer_game_(reference, column));
        value += integer_game_(reference, reference);
        value -= integer_game_(row, reference);
        reduced_entry_cache_ready_[offset] = 1;
    }
    return value;
}

/*
 * Let U=I\{m}; the indices in J\I are the outside best replies. With the unrestricted U coordinates first and the coordinates of
 * the outside best replies second, write the reduced extended Hessian as
 *
 *     R = [ H  G ],
 *         [ G' Q ]
 *
 * where dH, dG, and dQ are integer because integer_game_=dA. The retained factorization solves dH*N=dG*delta. Therefore
 *
 *     -(delta*dQ - dG'*N) = d*delta * [-(Q-G'*H^-1*G)].
 *
 * The bracketed matrix is one half of Bomze's final reduced matrix B^(r), so it differs only by an irrelevant positive scale. The
 * negation is retained deliberately: it lets every later check use the standard positive-definiteness and strict-copositivity signs
 * instead of reversing every comparison for a nonstandard "co-negativity" test. Original dG entries remain in
 * reduced_entry_cache_ after the solve overwrites the work matrix with N.
 */
size_t exact_candidate_solver::build_scaled_reduced_b(support::bitset support, support::bitset outside_best_replies,
                                                       numeric::matrix_int& result)
{
    assert(outside_best_replies != 0);

    uint8_t support_indices[support::kMaxBitsetDimension];
    uint8_t outside_indices[support::kMaxBitsetDimension];
    const size_t support_size = support::extract_set_indices(support, dimension_, support_indices);
    const size_t outside_size = support::extract_set_indices(outside_best_replies, dimension_, outside_indices);
    return build_scaled_reduced_b_from_indices(support_indices, support_size, outside_indices, outside_size, result);
}

size_t exact_candidate_solver::build_scaled_reduced_b(const support::bitset_multiword& support,
                                                   const support::bitset_multiword& outside_best_replies, numeric::matrix_int& result)
{
    assert(support.dimension() == dimension_);
    assert(outside_best_replies.dimension() == dimension_);
    assert(!outside_best_replies.empty());

    support.extract_set_indices(support_indices_large_);
    outside_best_replies.extract_set_indices(non_support_indices_large_);
    return build_scaled_reduced_b_from_indices(support_indices_large_.data(), support_indices_large_.size(),
                                               non_support_indices_large_.data(), non_support_indices_large_.size(), result);
}

template<class Index>
size_t exact_candidate_solver::build_scaled_reduced_b_from_indices(const Index* support_indices, size_t support_size,
                                                                const Index* outside_indices, size_t outside_size,
                                                                numeric::matrix_int& result)
{
    assert(reduced_hessian_is_negative_definite_);
    assert(outside_size > 0);

    const size_t reduced_dimension = support_size - 1;
    const size_t reference = support_indices[0];

    if (reduced_dimension > 0) {
        assert(reduced_system_dimension_ == reduced_dimension);
        stability_solution_numerators_.resize(reduced_dimension, outside_size);

        // Before the solve this matrix is dG. Afterwards it contains the solution numerators N.
        for (size_t row = 0; row < reduced_dimension; ++row) {
            const size_t i = support_indices[row + 1];
            for (size_t column = 0; column < outside_size; ++column) {
                stability_solution_numerators_(row, column) = reduced_entry(reference, i, outside_indices[column]);
            }
        }

        ffldlt_workspace_.solve_factored_negative_definite_inplace(stability_solution_numerators_, solution_denominator_, reduced_system_);
    }

    result.resize(outside_size, outside_size);

    numeric::integer scaled_entry;
    for (size_t row = 0; row < outside_size; ++row) {
        for (size_t column = 0; column <= row; ++column) {
            scaled_entry.set_product(
                solution_denominator_, reduced_entry(reference, outside_indices[row], outside_indices[column]));

            for (size_t inner = 0; inner < reduced_dimension; ++inner) {
                const size_t i = support_indices[inner + 1];
                scaled_entry.submul(
                    reduced_entry(reference, i, outside_indices[row]), stability_solution_numerators_(inner, column));
            }

            scaled_entry.negate();
            result(row, column) = scaled_entry;
            if (row != column) result(column, row) = result(row, column);
        }
    }
    return outside_size;
}

template<class Index>
void exact_candidate_solver::calculate_integer_payoff(numeric::integer& value, size_t strategy, size_t reference,
                                                   const Index* support_indices, size_t reduced_dimension)
{
    value.set_product(integer_game_(strategy, reference), reference_numerator_);
    for (size_t position = 0; position < reduced_dimension; ++position) {
        value.addmul(integer_game_(strategy, support_indices[position + 1]), solution_numerators_(position, 0));
    }
}

template<class SupportMask>
void exact_candidate_solver::ensure_candidate_vector(basic_candidate<SupportMask>& result) const
{
    if (result.vector.size() != dimension_) result.vector.resize(dimension_);
}

/*
 * A mixed strategy x with support S is a symmetric Nash equilibrium exactly when
 *
 *     A_S*x = u*1,          every used strategy earns the same payoff u,
 *     1^T*x = 1,            probabilities sum to one,
 *     x_i > 0 for i in S,   S is the actual support,
 *     (A*x)_i <= u outside S.
 *
 * The fraction-free LDL^T solve proves nonsingularity, returns all probabilities with one common denominator, and records the exact inertia of H.
 * A failed test returns immediately; rational candidate fields are materialized only after every exact candidate condition succeeds.
 */
bool exact_candidate_solver::find(const support::bitset& support, size_t support_size, candidate& result, bool materialize_output)
{
    uint8_t support_indices[support::kMaxBitsetDimension];
    uint8_t non_support_indices[support::kMaxBitsetDimension];
    const size_t support_count = support::extract_set_indices(support, dimension_, support_indices);
    const support::bitset complement = support::set_all_n_bits(dimension_) & ~support;
    const size_t non_support_count = support::extract_set_indices(complement, dimension_, non_support_indices);
    return find_from_indices(support, support_size, result, materialize_output, support_indices, support_count,
                             non_support_indices, non_support_count);
}

bool exact_candidate_solver::find(const support::bitset_multiword& support, size_t support_size, candidate_multiword& result,
                               bool materialize_output)
{
    assert(support.dimension() == dimension_);
    support.extract_set_indices(support_indices_large_);
    support.extract_unset_indices(non_support_indices_large_);

    return find_from_indices(support, support_size, result, materialize_output, support_indices_large_.data(),
                             support_indices_large_.size(), non_support_indices_large_.data(), non_support_indices_large_.size());
}

template<class SupportMask, class Index>
bool exact_candidate_solver::find_from_indices(const SupportMask& support, size_t support_size,
                                            basic_candidate<SupportMask>& result, bool materialize_output,
                                            const Index* support_indices, size_t support_count,
                                            const Index* non_support_indices, size_t non_support_count)
{
    reduced_hessian_is_negative_definite_ = false;
    reduced_hessian_negative_inertia_ = 0;

    assert(support_count == support_size);
    static_cast<void>(support_size); // Support generators guarantee this invariant when assertions are disabled.

    const size_t reference = support_indices[0];
    const size_t reduced_dimension = support_count - 1;

    if (reduced_dimension == 0) {
        // A pure support has no tangent direction. Its reduced Hessian is vacuously negative definite and normalization fixes its probability at one.
        solution_denominator_.set_one();
        reference_numerator_.set_one();
        reduced_hessian_is_negative_definite_ = true;
    } else {
        resize_reduced_system(reduced_dimension);
        build_reduced_system(support_indices, reduced_dimension);

        numeric::fraction_free_ldlt_inertia inertia;
        const bool nonsingular = ffldlt_workspace_.solve_inplace(
            solution_numerators_, solution_denominator_, reduced_system_, right_hand_side_, inertia);
        reduced_hessian_negative_inertia_ = static_cast<size_t>(inertia.negative);
        if (!nonsingular) return false;
        reduced_hessian_is_negative_definite_ = reduced_hessian_negative_inertia_ == reduced_dimension;

        // A3 produces ESS output only. Failed strict concavity is already a complete rejection and an exact upper-cone certificate,
        // so probabilities and outside payoffs cannot affect the result and need not be calculated.
        if (!reduced_hessian_is_negative_definite_) return false;

        // The solved entries are the probabilities for S without the reference. Recover the reference probability from x_m=1-sum(y).
        reference_numerator_ = solution_denominator_;
        for (size_t position = 0; position < reduced_dimension; ++position) {
            const auto numerator = solution_numerators_(position, 0);
            if (numerator.sign() <= 0) return false;
            reference_numerator_ -= numerator;
        }
        if (reference_numerator_.sign() <= 0) return false;
    }

    // Recover the common support payoff from the reference row. Its denominator is game_denominator*solution_denominator.
    calculate_integer_payoff(payoff_numerator_, reference, reference, support_indices, reduced_dimension);

    result.extended_support = support;
    for (size_t position = 0; position < non_support_count; ++position) {
        const size_t outside_strategy = non_support_indices[position];
        calculate_integer_payoff(outside_payoff_numerator_, outside_strategy, reference, support_indices, reduced_dimension);
        const int comparison = outside_payoff_numerator_.compare(payoff_numerator_);
        if (comparison > 0) return false;
        if (comparison == 0) {
            if constexpr (std::is_same_v<SupportMask, support::bitset>)
                result.extended_support = support::set_bit_at_pos(result.extended_support, outside_strategy);
            else
                result.extended_support.set_bit_at_pos(outside_strategy);
        }
    }

    if constexpr (std::is_same_v<SupportMask, support::bitset>)
        result.extended_support_size = support::count_set_bits(result.extended_support);
    else
        result.extended_support_size = result.extended_support.count_set_bits();

    // Stability, pruning, counting, and logging consume only exact integers and masks. Construct public rational output only when requested.
    if (materialize_output) {
        payoff_denominator_.set_product(game_denominator_, solution_denominator_);
        result.payoff.set_ratio(payoff_numerator_, payoff_denominator_);
        result.payoff_dbl = result.payoff.to_dbl();

        ensure_candidate_vector(result);
        for (size_t position = 0; position < non_support_count; ++position) result.vector[non_support_indices[position]].set_zero();
        result.vector[reference].set_ratio(reference_numerator_, solution_denominator_);
        for (size_t position = 0; position < reduced_dimension; ++position) {
            result.vector[support_indices[position + 1]].set_ratio(solution_numerators_(position, 0), solution_denominator_);
        }
    }

    return true;
}

} // namespace fracessa::search
