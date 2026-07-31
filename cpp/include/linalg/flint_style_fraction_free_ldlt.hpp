#pragma once

/*
 * The immediate-integer Bareiss update below is adapted from FLINT's fmpz_mat_fflu:
 *
 *     Copyright (C) 2011 Fredrik Johansson
 *
 * FLINT is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License (LGPL) as published by
 * the Free Software Foundation; either version 3 of the License, or (at your option) any later version. See <https://www.gnu.org/licenses/>.
 */

#include <flint/fmpz_mat.h>
#include <flint/ulong_extras.h>

#include <cassert>
#include <cstddef>
#include <vector>

namespace linalg {

struct fraction_free_ldlt_inertia {
    slong positive = 0;
    slong negative = 0;
};

/*
 * Reusable scratch for the in-place fraction-free symmetric solve. Keeping this object beside the caller's reusable matrices avoids allocations
 * inside a support-enumeration loop. One workspace must not be used concurrently by multiple solves.
 */
class fraction_free_ldlt_workspace {
public:
    explicit fraction_free_ldlt_workspace(size_t maximum_dimension)
        : coordinate_operations_(2 * maximum_dimension)
    {
        fmpz_init(previous_pivot_);
        fmpz_init(temporary_1_);
        fmpz_init(temporary_2_);
    }

    fraction_free_ldlt_workspace(const fraction_free_ldlt_workspace&) = delete;
    fraction_free_ldlt_workspace& operator=(const fraction_free_ldlt_workspace&) = delete;

    ~fraction_free_ldlt_workspace()
    {
        fmpz_clear(temporary_2_);
        fmpz_clear(temporary_1_);
        fmpz_clear(previous_pivot_);
    }

    /*
     * Solve system*X = right_hand_side*denominator for a nonsingular symmetric integer system while recording its exact inertia. The lower triangle
     * of system is read and overwritten by the fraction-free factorization; right_hand_side is also overwritten. X and right_hand_side currently
     * have one column. Returns 1 for a nonsingular system and 0 for a singular system. Outputs are undefined after a zero return.
     */
    int solve_inplace(fmpz_mat_t X, fmpz_t denominator, fmpz_mat_t system, fmpz_mat_t right_hand_side,
                      fraction_free_ldlt_inertia& inertia)
    {
        assert(fmpz_mat_nrows(system) == fmpz_mat_ncols(system));
        assert(fmpz_mat_nrows(right_hand_side) == fmpz_mat_nrows(system));
        assert(fmpz_mat_ncols(right_hand_side) == 1);
        assert(fmpz_mat_nrows(X) == fmpz_mat_nrows(system));
        assert(fmpz_mat_ncols(X) == 1);

        const size_t dimension = static_cast<size_t>(fmpz_mat_nrows(system));

        inertia = {};
        operation_count_ = 0;
        assert(coordinate_operations_.size() >= 2 * dimension);

        if (dimension == 0) {
            fmpz_one(denominator);
            return 1;
        }

        bool immediate = all_values_are_immediate(system, right_hand_side, dimension);
        bool previous_is_one = true;
        ulong unsigned_denominator = 0;
        ulong denominator_inverse = 0;
        slong normalization_shift = 0;
        bool denominator_is_negative = false;
        int previous_pivot_sign = 1;
        slong positive_inertia = 0;

        for (size_t pivot_position = 0; pivot_position < dimension; ++pivot_position) {
            if (!select_nonzero_diagonal(system, right_hand_side, pivot_position, dimension, immediate)) return 0;
            fmpz* pivot = entry(system, pivot_position, pivot_position);

            // Bareiss pivot p_k divided by p_(k-1) is the corresponding diagonal entry of ordinary LDL^T, so its sign gives one inertia count.
            const int diagonal_sign = fmpz_sgn(pivot) * previous_pivot_sign;
            if (diagonal_sign > 0) ++positive_inertia;
            previous_pivot_sign = fmpz_sgn(pivot);

            if (pivot_position + 1 == dimension) break;

            bool next_step_is_immediate = immediate;
            if (immediate) {
                // Every source operand was immediate when this pivot step began. A large result only changes subsequent steps, so this step may
                // continue using FLINT's two-limb update without rechecking every destination inside the loops.
                for (size_t row = pivot_position + 1; row < dimension; ++row) {
                    for (size_t column = pivot_position + 1; column <= row; ++column) {
                        next_step_is_immediate &= update_immediate(entry(system, row, column), pivot, entry(system, row, pivot_position),
                                                                   entry(system, column, pivot_position), pivot_position > 0, previous_is_one,
                                                                   unsigned_denominator, denominator_inverse, normalization_shift,
                                                                   denominator_is_negative);
                    }
                    next_step_is_immediate &= update_immediate(rhs_entry(right_hand_side, row), pivot, entry(system, row, pivot_position),
                                                               rhs_entry(right_hand_side, pivot_position), pivot_position > 0, previous_is_one,
                                                               unsigned_denominator, denominator_inverse, normalization_shift,
                                                               denominator_is_negative);
                }
            } else {
                for (size_t row = pivot_position + 1; row < dimension; ++row) {
                    for (size_t column = pivot_position + 1; column <= row; ++column) {
                        fmpz* result = entry(system, row, column);
                        fmpz_mul(result, result, pivot);
                        fmpz_submul(result, entry(system, row, pivot_position), entry(system, column, pivot_position));
                        if (pivot_position > 0 && !previous_is_one) fmpz_divexact(result, result, previous_pivot_);
                    }

                    fmpz* result = rhs_entry(right_hand_side, row);
                    fmpz_mul(result, result, pivot);
                    fmpz_submul(result, entry(system, row, pivot_position), rhs_entry(right_hand_side, pivot_position));
                    if (pivot_position > 0 && !previous_is_one) fmpz_divexact(result, result, previous_pivot_);
                }
            }

            fmpz_set(previous_pivot_, pivot);
            previous_is_one = fmpz_is_one(previous_pivot_);
            immediate = next_step_is_immediate;

            if (immediate) {
                unsigned_denominator = FLINT_ABS(static_cast<slong>(*previous_pivot_));
                denominator_is_negative = static_cast<slong>(*previous_pivot_) < 0;
                normalization_shift = flint_clz(unsigned_denominator);
                denominator_inverse = n_preinvert_limb_prenorm(unsigned_denominator << normalization_shift);
            }
        }

        inertia.positive = positive_inertia;
        inertia.negative = static_cast<slong>(dimension) - positive_inertia;
        fmpz_set(denominator, entry(system, dimension - 1, dimension - 1));

        // Back substitution stays integral: with det(system) as the common denominator, every solution numerator is a Cramer determinant.
        for (size_t row = dimension; row-- > 0; ) {
            fmpz* numerator = solution_entry(X, row);
            fmpz_mul(numerator, denominator, rhs_entry(right_hand_side, row));
            for (size_t column = row + 1; column < dimension; ++column) {
                fmpz_submul(numerator, entry(system, column, row), solution_entry(X, column));
            }
            fmpz_divexact(numerator, numerator, entry(system, row, row));
        }

        restore_original_coordinates(X);
        return 1;
    }

private:
    enum class coordinate_operation_kind : unsigned char { swap, add };

    struct coordinate_operation {
        coordinate_operation_kind kind;
        size_t target;
        size_t source;
    };

    static ulong right_shift(ulong value, unsigned int count) noexcept
    {
        return count == FLINT_BITS ? UWORD(0) : value >> count;
    }

    static fmpz* entry(fmpz_mat_struct* system, size_t row, size_t column) noexcept
    {
        return fmpz_mat_entry(system, static_cast<slong>(row), static_cast<slong>(column));
    }

    static fmpz* lower_entry(fmpz_mat_struct* system, size_t first, size_t second) noexcept
    {
        return first >= second ? entry(system, first, second) : entry(system, second, first);
    }

    static fmpz* rhs_entry(fmpz_mat_struct* right_hand_side, size_t row) noexcept
    {
        return fmpz_mat_entry(right_hand_side, static_cast<slong>(row), 0);
    }

    static fmpz* solution_entry(fmpz_mat_struct* X, size_t row) noexcept
    {
        return fmpz_mat_entry(X, static_cast<slong>(row), 0);
    }

    void remember_operation(coordinate_operation_kind kind, size_t target, size_t source)
    {
        assert(operation_count_ < coordinate_operations_.size());
        coordinate_operations_[operation_count_++] = {kind, target, source};
    }

    /*
     * Exchange two active coordinates. Completed rows of the fraction-free upper factor are stored transposed in the completed lower columns, so
     * they receive the corresponding column exchange separately. The active lower triangle and right-hand side receive the congruent exchange.
     */
    void swap_symmetric_coordinates(fmpz_mat_struct* system, fmpz_mat_struct* right_hand_side, size_t completed, size_t dimension,
                                    size_t first, size_t second)
    {
        if (first == second) return;

        for (size_t column = 0; column < completed; ++column) fmpz_swap(entry(system, first, column), entry(system, second, column));
        fmpz_swap(rhs_entry(right_hand_side, first), rhs_entry(right_hand_side, second));
        fmpz_swap(entry(system, first, first), entry(system, second, second));

        for (size_t index = completed; index < dimension; ++index) {
            if (index == first || index == second) continue;
            fmpz_swap(lower_entry(system, first, index), lower_entry(system, second, index));
        }

        remember_operation(coordinate_operation_kind::swap, first, second);
    }

    /*
     * Apply the unimodular change old_coordinate = T*new_coordinate with T=I+e_source*e_target^T. Equivalently, add active column source to column
     * target and active row source to row target. The new target diagonal is A_tt+2*A_ts+A_ss and is therefore nonzero when the active diagonal is
     * zero and A_ts is the selected nonzero off-diagonal entry.
     */
    void add_symmetric_coordinate(fmpz_mat_struct* system, fmpz_mat_struct* right_hand_side, size_t completed, size_t dimension,
                                  size_t target, size_t source)
    {
        for (size_t column = 0; column < completed; ++column) {
            fmpz_add(entry(system, target, column), entry(system, target, column), entry(system, source, column));
        }
        fmpz_add(rhs_entry(right_hand_side, target), rhs_entry(right_hand_side, target), rhs_entry(right_hand_side, source));

        fmpz_set(temporary_1_, entry(system, target, target));
        fmpz_set(temporary_2_, lower_entry(system, target, source));

        for (size_t index = completed; index < dimension; ++index) {
            if (index == target || index == source) continue;
            fmpz_add(lower_entry(system, target, index), lower_entry(system, target, index), lower_entry(system, source, index));
        }

        fmpz_add(lower_entry(system, target, source), temporary_2_, entry(system, source, source));
        fmpz_add(entry(system, target, target), temporary_1_, temporary_2_);
        fmpz_add(entry(system, target, target), entry(system, target, target), temporary_2_);
        fmpz_add(entry(system, target, target), entry(system, target, target), entry(system, source, source));

        remember_operation(coordinate_operation_kind::add, target, source);
    }

    static bool all_values_are_immediate(fmpz_mat_struct* system, fmpz_mat_struct* right_hand_side, size_t dimension) noexcept
    {
        for (size_t row = 0; row < dimension; ++row) {
            if (COEFF_IS_MPZ(*rhs_entry(right_hand_side, row))) return false;
            for (size_t column = 0; column <= row; ++column) {
                if (COEFF_IS_MPZ(*entry(system, row, column))) return false;
            }
        }
        return true;
    }

    /* Compute result=(result*pivot-left*right)/previous_pivot while every input is a FLINT immediate integer. */
    bool update_immediate(fmpz* result, const fmpz* pivot, const fmpz* left, const fmpz* right, bool divide_by_previous,
                          bool previous_is_one, ulong unsigned_denominator, ulong denominator_inverse, slong normalization_shift,
                          bool denominator_is_negative)
    {
        ulong first_high, first_low, second_high, second_low;
        smul_ppmm(first_high, first_low, *result, *pivot);
        smul_ppmm(second_high, second_low, *left, *right);
        sub_ddmmss(first_high, first_low, first_high, first_low, second_high, second_low);

        const bool result_is_negative = static_cast<slong>(first_high) < 0;
        if (result_is_negative) sub_ddmmss(first_high, first_low, UWORD(0), UWORD(0), first_high, first_low);

        if (divide_by_previous && !previous_is_one) {
            if (first_high >= unsigned_denominator) {
                fmpz_set_uiui(result, first_high, first_low);
                if (result_is_negative) fmpz_neg(result, result);
                fmpz_divexact(result, result, previous_pivot_);
            } else {
                ulong quotient;
                ulong FLINT_SET_BUT_UNUSED(remainder);
                udiv_qrnnd_preinv(quotient, remainder,
                                  (first_high << normalization_shift) + right_shift(first_low, FLINT_BITS - normalization_shift),
                                  first_low << normalization_shift, unsigned_denominator << normalization_shift, denominator_inverse);
                if (result_is_negative != denominator_is_negative) fmpz_neg_ui(result, quotient);
                else fmpz_set_ui(result, quotient);
            }
        } else {
            if (first_high > 0) fmpz_set_uiui(result, first_high, first_low);
            else fmpz_set_ui(result, first_low);
            if (result_is_negative) fmpz_neg(result, result);
        }

        return !COEFF_IS_MPZ(*result);
    }

    bool select_nonzero_diagonal(fmpz_mat_struct* system, fmpz_mat_struct* right_hand_side, size_t pivot_position, size_t dimension,
                                 bool& immediate)
    {
        size_t diagonal_pivot = pivot_position;
        while (diagonal_pivot < dimension && fmpz_is_zero(entry(system, diagonal_pivot, diagonal_pivot))) ++diagonal_pivot;

        if (diagonal_pivot < dimension) {
            swap_symmetric_coordinates(system, right_hand_side, pivot_position, dimension, pivot_position, diagonal_pivot);
            return true;
        }

        size_t first = dimension;
        size_t second = dimension;
        for (size_t row = pivot_position + 1; row < dimension && first == dimension; ++row) {
            for (size_t column = pivot_position; column < row; ++column) {
                if (!fmpz_is_zero(entry(system, row, column))) {
                    first = column;
                    second = row;
                    break;
                }
            }
        }
        if (first == dimension) return false;

        swap_symmetric_coordinates(system, right_hand_side, pivot_position, dimension, pivot_position, first);
        add_symmetric_coordinate(system, right_hand_side, pivot_position, dimension, pivot_position, second);
        if (immediate) immediate = all_values_are_immediate(system, right_hand_side, dimension);
        assert(!fmpz_is_zero(entry(system, pivot_position, pivot_position)));
        return true;
    }

    void restore_original_coordinates(fmpz_mat_struct* X)
    {
        for (size_t index = operation_count_; index-- > 0; ) {
            const coordinate_operation operation = coordinate_operations_[index];
            if (operation.kind == coordinate_operation_kind::swap) {
                fmpz_swap(solution_entry(X, operation.target), solution_entry(X, operation.source));
            } else {
                // The forward congruence used old=T*new with T=I+e_source*e_target^T, hence x_source <- x_source+x_target.
                fmpz_add(solution_entry(X, operation.source), solution_entry(X, operation.source), solution_entry(X, operation.target));
            }
        }
    }

    size_t operation_count_ = 0;
    std::vector<coordinate_operation> coordinate_operations_;
    fmpz_t previous_pivot_;
    fmpz_t temporary_1_;
    fmpz_t temporary_2_;
};

/* Fast form for callers that already own disposable matrices and reusable scratch. */
inline int fmpz_mat_solve_ffldlt_inplace(fmpz_mat_t X, fmpz_t denominator, fmpz_mat_t system, fmpz_mat_t right_hand_side,
                                         fraction_free_ldlt_inertia& inertia, fraction_free_ldlt_workspace& workspace)
{
    return workspace.solve_inplace(X, denominator, system, right_hand_side, inertia);
}

/*
 * FLINT-style convenience interface. Like fmpz_mat_solve_fflu, this preserves the const inputs and computes integer X and denominator satisfying
 * system*X = right_hand_side*denominator. Repeated hot-path callers should use the in-place overload above to avoid copying and allocation.
 */
inline int fmpz_mat_solve_ffldlt(fmpz_mat_t X, fmpz_t denominator, const fmpz_mat_t system, const fmpz_mat_t right_hand_side,
                                 fraction_free_ldlt_inertia& inertia)
{
    const size_t dimension = static_cast<size_t>(fmpz_mat_nrows(system));
    fraction_free_ldlt_workspace workspace(dimension);
    fmpz_mat_t system_copy;
    fmpz_mat_t right_hand_side_copy;
    fmpz_mat_init_set(system_copy, system);
    fmpz_mat_init_set(right_hand_side_copy, right_hand_side);

    const int nonsingular = workspace.solve_inplace(X, denominator, system_copy, right_hand_side_copy, inertia);
    fmpz_mat_clear(right_hand_side_copy);
    fmpz_mat_clear(system_copy);
    return nonsingular;
}

} // namespace linalg
