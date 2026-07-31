#pragma once

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>

#include <fracessa/bitset64.hpp>
#include <fracessa/candidate.hpp>
#include <linalg/flint_style_fraction_free_ldlt.hpp>
#include <linalg/matrix_fraction.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>

/*
 * Exact candidate solver specialized for FracESSA's border-reduced symmetric system.
 *
 * For support S={s_0,...,s_k}, write x=e_(s_0)+Z*y with columns e_(s_i)-e_(s_0), i>0. The bordered candidate equations reduce to
 *
 *     H*y = r,
 *     H_ij = A_(s_i,s_j)-A_(s_i,s_0)-A_(s_0,s_j)+A_(s_0,s_0),
 *     r_i  = A_(s_0,s_0)-A_(s_i,s_0).
 *
 * The game is converted once to B=d*A with one positive common denominator d. Hence dH and dr are integral, have the same solution as H*y=r,
 * and dH has the same inertia as H. Symmetric Bareiss elimination then performs one fraction-free LDL^T-style factorization which provides all
 * three facts needed by candidate search:
 *
 *   - a missing pivot proves singularity;
 *   - fraction-free back substitution returns the exact candidate with one common denominator;
 *   - the signs of consecutive Bareiss pivots give the inertia of H.
 *
 * Exact arithmetic needs no numerical pivot threshold. A nonzero active diagonal is a valid 1x1 pivot. If an active nonsingular symmetric block
 * has zero diagonal, it must have a nonzero off-diagonal entry a. The unimodular coordinate change e_i <- e_i+e_j changes that diagonal to 2a.
 * Applying the same change as a row and a column operation preserves symmetry, determinant, inertia, and exact divisibility while avoiding a
 * separate fraction-free 2x2-block implementation. The coordinate changes are undone after back substitution.
 */
class ffldlt_candidate_solver {
public:
    explicit ffldlt_candidate_solver(const linalg::matrix_frc& game)
        : dimension_(game.rows()), ffldlt_workspace_(dimension_)
    {
        fmpq_mat_t rational_game;
        fmpq_mat_init(rational_game, static_cast<slong>(dimension_), static_cast<slong>(dimension_));
        fmpz_mat_init(integer_game_, static_cast<slong>(dimension_), static_cast<slong>(dimension_));
        fmpz_init(game_denominator_);

        for (size_t row = 0; row < dimension_; ++row) {
            for (size_t column = 0; column < dimension_; ++column) {
                fmpq_set(fmpq_mat_entry(rational_game, static_cast<slong>(row), static_cast<slong>(column)), raw_fraction(game(row, column)));
            }
        }
        fmpq_mat_get_fmpz_mat_matwise(integer_game_, game_denominator_, rational_game);
        fmpq_mat_clear(rational_game);

        fmpz_mat_init(system_, 0, 0);
        fmpz_mat_init(right_hand_side_, 0, 0);
        fmpz_mat_init(solution_numerators_, 0, 0);
        fmpz_init(solution_denominator_);
        fmpz_init(reference_numerator_);
        fmpz_init(payoff_numerator_);
        fmpz_init(payoff_denominator_);
        fmpz_init(outside_sum_);
    }

    ffldlt_candidate_solver(const ffldlt_candidate_solver&) = delete;
    ffldlt_candidate_solver& operator=(const ffldlt_candidate_solver&) = delete;

    ~ffldlt_candidate_solver()
    {
        fmpz_clear(outside_sum_);
        fmpz_clear(payoff_denominator_);
        fmpz_clear(payoff_numerator_);
        fmpz_clear(reference_numerator_);
        fmpz_clear(solution_denominator_);
        fmpz_mat_clear(solution_numerators_);
        fmpz_mat_clear(right_hand_side_);
        fmpz_mat_clear(system_);
        fmpz_clear(game_denominator_);
        fmpz_mat_clear(integer_game_);
    }

    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
    {
        reduced_hessian_is_negative_definite_ = false;

        uint8_t support_indices[bs64::kMaxBitsetDimension];
        uint8_t non_support_indices[bs64::kMaxBitsetDimension];
        const size_t support_count = bs64::extract_set_indices(support, dimension_, support_indices);
        const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
        const size_t non_support_count = bs64::extract_set_indices(complement, dimension_, non_support_indices);
        assert(support_count == support_size);
        static_cast<void>(support_size); // The assertion is compiled out in Release; support generators guarantee this invariant.

        const size_t reference = support_indices[0];
        const size_t reduced_dimension = support_count - 1;

        if (reduced_dimension == 0) {
            fmpz_one(solution_denominator_);
            fmpz_one(reference_numerator_);
            reduced_hessian_is_negative_definite_ = true;
        } else {
            resize_reduced_system(reduced_dimension);
            build_reduced_system(support_indices, reduced_dimension);
            linalg::fraction_free_ldlt_inertia inertia;
            if (!linalg::fmpz_mat_solve_ffldlt_inplace(solution_numerators_, solution_denominator_, system_, right_hand_side_, inertia,
                                                        ffldlt_workspace_)) return false;
            reduced_hessian_is_negative_definite_ = inertia.positive == 0;

            if (fmpz_sgn(solution_denominator_) < 0) {
                fmpz_neg(solution_denominator_, solution_denominator_);
                fmpz_mat_neg(solution_numerators_, solution_numerators_);
            }

            fmpz_set(reference_numerator_, solution_denominator_);
            for (size_t position = 0; position < reduced_dimension; ++position) {
                const fmpz* numerator = solution_entry(position);
                if (fmpz_sgn(numerator) <= 0) return false;
                fmpz_sub(reference_numerator_, reference_numerator_, numerator);
            }
            if (fmpz_sgn(reference_numerator_) <= 0) return false;
        }

        calculate_payoff(reference, support_indices, reduced_dimension);

        result.extended_support = support;
        for (size_t position = 0; position < non_support_count; ++position) {
            const size_t outside_strategy = non_support_indices[position];
            calculate_integer_payoff(outside_sum_, outside_strategy, reference, support_indices, reduced_dimension);
            const int comparison = fmpz_cmp(outside_sum_, payoff_numerator_);
            if (comparison > 0) return false;
            if (comparison == 0) result.extended_support = bs64::set_bit_at_pos(result.extended_support, outside_strategy);
        }

        fmpz_mul(payoff_denominator_, game_denominator_, solution_denominator_);
        fmpq_set_fmpz_frac(result.payoff.data(), payoff_numerator_, payoff_denominator_);
        result.payoff_dbl = result.payoff.to_dbl();
        result.extended_support_size = bs64::count_set_bits(result.extended_support);

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

    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

private:
    static const fmpq* raw_fraction(const fraction& value) noexcept
    {
        return const_cast<fraction&>(value).data();
    }

    fmpz* entry(size_t row, size_t column) noexcept
    {
        return fmpz_mat_entry(system_, static_cast<slong>(row), static_cast<slong>(column));
    }

    fmpz* rhs_entry(size_t row) noexcept
    {
        return fmpz_mat_entry(right_hand_side_, static_cast<slong>(row), 0);
    }

    fmpz* solution_entry(size_t row) noexcept
    {
        return fmpz_mat_entry(solution_numerators_, static_cast<slong>(row), 0);
    }

    const fmpz* solution_entry(size_t row) const noexcept
    {
        return fmpz_mat_entry(solution_numerators_, static_cast<slong>(row), 0);
    }

    const fmpz* game_entry(size_t row, size_t column) const noexcept
    {
        return fmpz_mat_entry(integer_game_, static_cast<slong>(row), static_cast<slong>(column));
    }

    void resize_reduced_system(size_t reduced_dimension)
    {
        if (reduced_dimension_ == reduced_dimension) return;

        fmpz_mat_clear(solution_numerators_);
        fmpz_mat_clear(right_hand_side_);
        fmpz_mat_clear(system_);

        reduced_dimension_ = reduced_dimension;
        const slong size = static_cast<slong>(reduced_dimension_);
        fmpz_mat_init(system_, size, size);
        fmpz_mat_init(right_hand_side_, size, 1);
        fmpz_mat_init(solution_numerators_, size, 1);
    }

    void build_reduced_system(const uint8_t* support_indices, size_t reduced_dimension)
    {
        const size_t reference = support_indices[0];
        const fmpz* reference_diagonal = game_entry(reference, reference);

        for (size_t row = 0; row < reduced_dimension; ++row) {
            const size_t i = support_indices[row + 1];
            fmpz_sub(rhs_entry(row), reference_diagonal, game_entry(i, reference));

            for (size_t column = 0; column <= row; ++column) {
                const size_t j = support_indices[column + 1];
                fmpz* value = entry(row, column);
                fmpz_set(value, game_entry(i, j));
                fmpz_sub(value, value, game_entry(i, reference));
                fmpz_sub(value, value, game_entry(reference, j));
                fmpz_add(value, value, reference_diagonal);
            }
        }
    }

    void calculate_integer_payoff(fmpz* value, size_t strategy, size_t reference, const uint8_t* support_indices, size_t reduced_dimension)
    {
        fmpz_mul(value, game_entry(strategy, reference), reference_numerator_);
        for (size_t position = 0; position < reduced_dimension; ++position) {
            fmpz_addmul(value, game_entry(strategy, support_indices[position + 1]), solution_entry(position));
        }
    }

    void calculate_payoff(size_t reference, const uint8_t* support_indices, size_t reduced_dimension)
    {
        calculate_integer_payoff(payoff_numerator_, reference, reference, support_indices, reduced_dimension);
    }

    void ensure_candidate_vector(candidate& result) const
    {
        if (result.vector.rows() != dimension_ || result.vector.cols() != 1) result.vector = linalg::matrix_frc(dimension_, 1);
    }

    size_t dimension_;
    size_t reduced_dimension_ = 0;
    bool reduced_hessian_is_negative_definite_ = false;
    linalg::fraction_free_ldlt_workspace ffldlt_workspace_;
    fmpz_mat_t integer_game_;
    fmpz_t game_denominator_;
    fmpz_mat_t system_;
    fmpz_mat_t right_hand_side_;
    fmpz_mat_t solution_numerators_;
    fmpz_t solution_denominator_;
    fmpz_t reference_numerator_;
    fmpz_t payoff_numerator_;
    fmpz_t payoff_denominator_;
    fmpz_t outside_sum_;
};
