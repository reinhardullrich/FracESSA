#include <fracessa/find_candidate_exact.hpp>

#include <flint/fmpq_mat.h>

#include <cassert>
#include <cstdint>

namespace candidate_search {
namespace {

const fmpq* raw_fraction(const fraction& value) noexcept
{
    return const_cast<fraction&>(value).data();
}

} // namespace

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
find_candidate_exact::find_candidate_exact(const linalg::matrix_frc& game_matrix)
    : dimension_(game_matrix.rows())
    , ffldlt_workspace_(dimension_)
{
    fmpq_mat_t rational_game;
    fmpq_mat_init(rational_game, static_cast<slong>(dimension_), static_cast<slong>(dimension_));
    fmpz_mat_init(integer_game_, static_cast<slong>(dimension_), static_cast<slong>(dimension_));
    fmpz_init(game_denominator_);

    for (size_t row = 0; row < dimension_; ++row) {
        for (size_t column = 0; column < dimension_; ++column) {
            fmpq_set(fmpq_mat_entry(rational_game, static_cast<slong>(row), static_cast<slong>(column)), raw_fraction(game_matrix(row, column)));
        }
    }
    fmpq_mat_get_fmpz_mat_matwise(integer_game_, game_denominator_, rational_game);
    fmpq_mat_clear(rational_game);

    fmpz_mat_init(reduced_system_, 0, 0);
    fmpz_mat_init(right_hand_side_, 0, 0);
    fmpz_mat_init(solution_numerators_, 0, 0);
    fmpz_init(solution_denominator_);
    fmpz_init(reference_numerator_);
    fmpz_init(payoff_numerator_);
    fmpz_init(payoff_denominator_);
    fmpz_init(outside_payoff_numerator_);
}

find_candidate_exact::~find_candidate_exact()
{
    fmpz_clear(outside_payoff_numerator_);
    fmpz_clear(payoff_denominator_);
    fmpz_clear(payoff_numerator_);
    fmpz_clear(reference_numerator_);
    fmpz_clear(solution_denominator_);
    fmpz_mat_clear(solution_numerators_);
    fmpz_mat_clear(right_hand_side_);
    fmpz_mat_clear(reduced_system_);
    fmpz_clear(game_denominator_);
    fmpz_mat_clear(integer_game_);
}

fmpz* find_candidate_exact::reduced_entry(size_t row, size_t column) noexcept
{
    return fmpz_mat_entry(reduced_system_, static_cast<slong>(row), static_cast<slong>(column));
}

fmpz* find_candidate_exact::right_hand_side_entry(size_t row) noexcept
{
    return fmpz_mat_entry(right_hand_side_, static_cast<slong>(row), 0);
}

fmpz* find_candidate_exact::solution_entry(size_t row) noexcept
{
    return fmpz_mat_entry(solution_numerators_, static_cast<slong>(row), 0);
}

const fmpz* find_candidate_exact::solution_entry(size_t row) const noexcept
{
    return fmpz_mat_entry(solution_numerators_, static_cast<slong>(row), 0);
}

const fmpz* find_candidate_exact::game_entry(size_t row, size_t column) const noexcept
{
    return fmpz_mat_entry(integer_game_, static_cast<slong>(row), static_cast<slong>(column));
}

void find_candidate_exact::resize_reduced_system(size_t reduced_dimension)
{
    if (reduced_dimension_ == reduced_dimension) return;

    fmpz_mat_clear(solution_numerators_);
    fmpz_mat_clear(right_hand_side_);
    fmpz_mat_clear(reduced_system_);

    reduced_dimension_ = reduced_dimension;
    const slong size = static_cast<slong>(reduced_dimension_);
    fmpz_mat_init(reduced_system_, size, size);
    fmpz_mat_init(right_hand_side_, size, 1);
    fmpz_mat_init(solution_numerators_, size, 1);
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
void find_candidate_exact::build_reduced_system(const uint8_t* support_indices, size_t reduced_dimension)
{
    const size_t reference = support_indices[0];
    const fmpz* reference_diagonal = game_entry(reference, reference);

    for (size_t row = 0; row < reduced_dimension; ++row) {
        const size_t i = support_indices[row + 1];
        fmpz_sub(right_hand_side_entry(row), reference_diagonal, game_entry(i, reference));

        // The symmetric fraction-free factorization reads only the lower triangle.
        for (size_t column = 0; column <= row; ++column) {
            const size_t j = support_indices[column + 1];
            fmpz* value = reduced_entry(row, column);
            fmpz_set(value, game_entry(i, j));
            fmpz_sub(value, value, game_entry(i, reference));
            fmpz_sub(value, value, game_entry(reference, j));
            fmpz_add(value, value, reference_diagonal);
        }
    }
}

void find_candidate_exact::calculate_integer_payoff(fmpz* value, size_t strategy, size_t reference, const uint8_t* support_indices,
                                                     size_t reduced_dimension)
{
    fmpz_mul(value, game_entry(strategy, reference), reference_numerator_);
    for (size_t position = 0; position < reduced_dimension; ++position) {
        fmpz_addmul(value, game_entry(strategy, support_indices[position + 1]), solution_entry(position));
    }
}

void find_candidate_exact::ensure_candidate_vector(candidate& result) const
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
bool find_candidate_exact::find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
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
        fmpz_one(solution_denominator_);
        fmpz_one(reference_numerator_);
        reduced_hessian_is_negative_definite_ = true;
    } else {
        resize_reduced_system(reduced_dimension);
        build_reduced_system(support_indices, reduced_dimension);

        linalg::fraction_free_ldlt_inertia inertia;
        if (!linalg::fmpz_mat_solve_ffldlt_inplace(solution_numerators_, solution_denominator_, reduced_system_, right_hand_side_, inertia,
                                                    ffldlt_workspace_)) return false;
        reduced_hessian_is_negative_definite_ = inertia.positive == 0;

        // FLINT rationals require a positive denominator. Negating both sides leaves every exact solution value unchanged.
        if (fmpz_sgn(solution_denominator_) < 0) {
            fmpz_neg(solution_denominator_, solution_denominator_);
            fmpz_mat_neg(solution_numerators_, solution_numerators_);
        }

        // The solved entries are the probabilities for S without the reference. Recover the reference probability from x_m=1-sum(y).
        fmpz_set(reference_numerator_, solution_denominator_);
        for (size_t position = 0; position < reduced_dimension; ++position) {
            const fmpz* numerator = solution_entry(position);
            if (fmpz_sgn(numerator) <= 0) return false;
            fmpz_sub(reference_numerator_, reference_numerator_, numerator);
        }
        if (fmpz_sgn(reference_numerator_) <= 0) return false;
    }

    // Recover the common support payoff from the reference row. Its denominator is game_denominator*solution_denominator.
    calculate_integer_payoff(payoff_numerator_, reference, reference, support_indices, reduced_dimension);

    result.extended_support = support;
    for (size_t position = 0; position < non_support_count; ++position) {
        const size_t outside_strategy = non_support_indices[position];
        calculate_integer_payoff(outside_payoff_numerator_, outside_strategy, reference, support_indices, reduced_dimension);
        const int comparison = fmpz_cmp(outside_payoff_numerator_, payoff_numerator_);
        if (comparison > 0) return false;
        if (comparison == 0) result.extended_support = bs64::set_bit_at_pos(result.extended_support, outside_strategy);
    }

    fmpz_mul(payoff_denominator_, game_denominator_, solution_denominator_);
    fmpq_set_fmpz_frac(result.payoff.data(), payoff_numerator_, payoff_denominator_);
    result.payoff_dbl = result.payoff.to_dbl();
    result.extended_support_size = bs64::count_set_bits(result.extended_support);

    // Stability does not consume the dense probability vector. Build it only when candidate output or logging requested it.
    if (materialize_vector) {
        ensure_candidate_vector(result);
        for (size_t position = 0; position < non_support_count; ++position) result.vector(non_support_indices[position], 0) = fraction::zero();
        fmpq_set_fmpz_frac(result.vector(reference, 0).data(), reference_numerator_, solution_denominator_);
        for (size_t position = 0; position < reduced_dimension; ++position) {
            fmpq_set_fmpz_frac(result.vector(support_indices[position + 1], 0).data(), solution_entry(position), solution_denominator_);
        }
    }

    return true;
}

} // namespace candidate_search
