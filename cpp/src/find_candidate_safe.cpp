#include <fracessa/find_candidate_safe.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

namespace candidate_search {

/*
 * Exact candidate solver for a symmetric payoff matrix.
 *
 * The constructor clears the rational denominators once for the whole game:
 *
 *     integer_game = game_denominator * game.
 *
 * The denominator is positive, so this scaling preserves every equality, inequality, and inertia sign used by candidate search and stability.
 * Each support can therefore stay in integer arithmetic until a successful candidate is written to the public rational result.
 */
find_candidate_safe::find_candidate_safe(const linalg::matrix_frc& game_matrix)
    : dimension_(game_matrix.rows())
    , ffldlt_workspace_(dimension_)
{
    integer_game_.set_from_fraction_matrix(game_matrix, game_denominator_);
}

/*
 * Let d be the least common positive denominator of the game and Z=d*A its integer matrix. The exact precision span is
 *
 *     P = M/m,
 *     M = max(d, max |Z_ij|),
 *     m = min(d, nonzero |Z_ij|, nonzero |Z_ij-Z_kl|).
 *
 * The scale d represents the 1 and -1 entries in the bordered double system. After sorting the symmetric upper triangle, only
 * adjacent distinct integers need to be compared to find the smallest nonzero pairwise difference.
 */
bool find_candidate_safe::precision_span_at_least(unsigned long limit) const
{
    const size_t entry_count = dimension_ * (dimension_ + 1) / 2;
    std::vector<linalg::integer::const_reference> entries;
    entries.reserve(entry_count);

    for (size_t row = 0; row < dimension_; ++row) {
        for (size_t column = row; column < dimension_; ++column) {
            entries.push_back(integer_game_(row, column));
        }
    }
    std::sort(entries.begin(), entries.end(), [](auto left, auto right) { return left.compare(right) < 0; });

    linalg::integer maximum(game_denominator_);
    linalg::integer minimum(game_denominator_);
    linalg::integer difference;
    linalg::integer scaled_minimum;

    for (const auto entry : entries) {
        if (entry.is_zero()) continue;
        if (entry.compare_abs(maximum) > 0) maximum.set_abs(entry);
        if (entry.compare_abs(minimum) < 0) minimum.set_abs(entry);
    }
    for (size_t index = 1; index < entries.size(); ++index) {
        if (entries[index - 1].compare(entries[index]) == 0) continue;
        difference.set_difference(entries[index], entries[index - 1]);
        if (difference.compare(minimum) < 0) minimum = difference;
    }

    scaled_minimum = minimum;
    scaled_minimum.multiply(limit);
    return maximum.compare(scaled_minimum) >= 0;
}

void find_candidate_safe::resize_reduced_system(size_t reduced_dimension)
{
    if (reduced_dimension_ == reduced_dimension) return;

    reduced_dimension_ = reduced_dimension;
    reduced_system_.resize(reduced_dimension_, reduced_dimension_);
    right_hand_side_.resize(reduced_dimension_, 1);
    solution_numerators_.resize(reduced_dimension_, 1);
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
void find_candidate_safe::build_reduced_system(const uint8_t* support_indices, size_t reduced_dimension)
{
    const size_t reference = support_indices[0];
    const auto reference_diagonal = integer_game_(reference, reference);

    for (size_t row = 0; row < reduced_dimension; ++row) {
        const size_t i = support_indices[row + 1];
        right_hand_side_(row, 0).set_difference(reference_diagonal, integer_game_(i, reference));

        // The symmetric fraction-free factorization reads only the lower triangle.
        for (size_t column = 0; column <= row; ++column) {
            const size_t j = support_indices[column + 1];
            auto value = reduced_system_(row, column);
            value = integer_game_(i, j);
            value -= integer_game_(i, reference);
            value -= integer_game_(reference, j);
            value += reference_diagonal;
        }
    }
}

void find_candidate_safe::calculate_integer_payoff(linalg::integer& value, size_t strategy, size_t reference,
                                                     const uint8_t* support_indices, size_t reduced_dimension)
{
    value.set_product(integer_game_(strategy, reference), reference_numerator_);
    for (size_t position = 0; position < reduced_dimension; ++position) {
        value.addmul(integer_game_(strategy, support_indices[position + 1]), solution_numerators_(position, 0));
    }
}

void find_candidate_safe::ensure_candidate_vector(candidate& result) const
{
    if (result.vector.rows() != dimension_ || result.vector.cols() != 1) result.vector = linalg::matrix_frc(dimension_, 1);
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
bool find_candidate_safe::find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
{
    reduced_hessian_is_negative_definite_ = false;

    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    const size_t support_count = bs64::extract_set_indices(support, dimension_, support_indices);
    const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
    const size_t non_support_count = bs64::extract_set_indices(complement, dimension_, non_support_indices);
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

        linalg::fraction_free_ldlt_inertia inertia;
        if (!ffldlt_workspace_.solve_inplace(solution_numerators_, solution_denominator_, reduced_system_, right_hand_side_, inertia)) return false;
        reduced_hessian_is_negative_definite_ = inertia.positive == 0;

        // FLINT rationals require a positive denominator. Negating both sides leaves every exact solution value unchanged.
        if (solution_denominator_.sign() < 0) {
            solution_denominator_.negate();
            solution_numerators_.negate();
        }

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
        if (comparison == 0) result.extended_support = bs64::set_bit_at_pos(result.extended_support, outside_strategy);
    }

    payoff_denominator_.set_product(game_denominator_, solution_denominator_);
    result.payoff.set_ratio(payoff_numerator_, payoff_denominator_);
    result.payoff_dbl = result.payoff.to_dbl();
    result.extended_support_size = bs64::count_set_bits(result.extended_support);

    // Stability does not consume the dense probability vector. Build it only when candidate output or logging requested it.
    if (materialize_vector) {
        ensure_candidate_vector(result);
        for (size_t position = 0; position < non_support_count; ++position) result.vector(non_support_indices[position], 0) = fraction::zero();
        result.vector(reference, 0).set_ratio(reference_numerator_, solution_denominator_);
        for (size_t position = 0; position < reduced_dimension; ++position) {
            result.vector(support_indices[position + 1], 0).set_ratio(solution_numerators_(position, 0), solution_denominator_);
        }
    }

    return true;
}

} // namespace candidate_search
