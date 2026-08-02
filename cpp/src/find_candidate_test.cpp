#include <fracessa/find_candidate_test.hpp>
#include <fracessa/find_candidate_safe.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

namespace candidate_search {
namespace {

constexpr double kPivotCutoff = 1e-12;
constexpr unsigned long kPrecisionSpanCutoff = 1'000'000'000UL;
constexpr double kProbabilityTolerance = -1e-10;
constexpr double kBunchKaufmanAlpha = 0.6403882032022076; // (1 + sqrt(17)) / 8
constexpr size_t kEquilibrationIterations = 100;

/*
 * The symmetric equilibration, indefinite factorization, and solve below adapt the lower-triangle paths of LAPACK 3.12.1
 * DSYEQUB, DSYTF2, and DSYTRS to FracESSA's small row-major matrices and single right-hand side. Changes from LAPACK are
 * zero-based indexing, direct loops in place of BLAS calls, fixed-size caller-owned arrays, and an inconclusive return for failed
 * equilibration or a small or non-finite pivot. The routine names and interfaces are local to this implementation, as requested by
 * LAPACK for modified routines.
 *
 * Copyright (c) 1992-2023 The University of Tennessee and The University of Tennessee Research Foundation. All rights reserved.
 * Copyright (c) 2000-2023 The University of California Berkeley. All rights reserved.
 * Copyright (c) 2006-2023 The University of Colorado Denver. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
 * conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided with the distribution.
 * - Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products
 *   derived from this software without specific prior written permission.
 *
 * The copyright holders provide no reassurances that the source code provided does not infringe any patent, copyright, or other
 * intellectual property rights of third parties. The copyright holders disclaim liability to any recipient for claims brought
 * by any third party for infringement of that party's intellectual property rights.
 *
 * THE SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 * BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Compute LAPACK's symmetric BIN equilibration once for the complete game and replace A by D*A*D. Every support then takes its
 * principal submatrix from the same scaled game. The transformed normalization vector D*1 is handled explicitly when the support
 * border is eliminated below, so this is a change of variables rather than a change of game.
 */
safe_fallback equilibrate_game_matrix(linalg::matrix_dbl& game, size_t dimension, double* scales)
{
    double beta[bs64::kMaxBitsetDimension];
    size_t active_coordinates[bs64::kMaxBitsetDimension];
    size_t active_dimension = 0;

    for (size_t coordinate = 0; coordinate < dimension; ++coordinate) scales[coordinate] = 0.0;
    for (size_t column = 0; column < dimension; ++column) {
        const double diagonal = std::abs(game(column, column));
        if (!std::isfinite(diagonal)) return safe_fallback::equilibration_invalid;
        scales[column] = std::max(scales[column], diagonal);
        for (size_t row = column + 1; row < dimension; ++row) {
            const double magnitude = std::abs(game(row, column));
            if (!std::isfinite(magnitude)) return safe_fallback::equilibration_invalid;
            scales[row] = std::max(scales[row], magnitude);
            scales[column] = std::max(scales[column], magnitude);
        }
    }

    for (size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        if (!std::isfinite(scales[coordinate])) return safe_fallback::equilibration_invalid;
        if (scales[coordinate] == 0.0) {
            // An exact zero row is also a zero column because the game is symmetric. Its positive scale is arbitrary, so leave it
            // unchanged and equilibrate every nonzero coordinate independently instead of rejecting or accepting the whole game.
            scales[coordinate] = 1.0;
            continue;
        }
        scales[coordinate] = 1.0 / scales[coordinate];
        if (!(scales[coordinate] > 0.0) || !std::isfinite(scales[coordinate])) return safe_fallback::equilibration_invalid;
        active_coordinates[active_dimension++] = coordinate;
    }
    if (active_dimension == 0) return safe_fallback::none;

    const double active_dimension_dbl = static_cast<double>(active_dimension);
    const double tolerance = 1.0 / std::sqrt(2.0 * active_dimension_dbl);
    double average = 0.0;
    bool converged = false;
    for (size_t iteration = 0; iteration < kEquilibrationIterations; ++iteration) {
        for (size_t active = 0; active < active_dimension; ++active) beta[active_coordinates[active]] = 0.0;
        for (size_t active_column = 0; active_column < active_dimension; ++active_column) {
            const size_t column = active_coordinates[active_column];
            beta[column] += std::abs(game(column, column)) * scales[column];
            for (size_t active_row = active_column + 1; active_row < active_dimension; ++active_row) {
                const size_t row = active_coordinates[active_row];
                const double magnitude = std::abs(game(row, column));
                beta[row] += magnitude * scales[column];
                beta[column] += magnitude * scales[row];
            }
        }

        average = 0.0;
        for (size_t active = 0; active < active_dimension; ++active) {
            const size_t coordinate = active_coordinates[active];
            average += scales[coordinate] * beta[coordinate];
        }
        average /= active_dimension_dbl;
        if (!(average > 0.0) || !std::isfinite(average)) return safe_fallback::equilibration_invalid;

        // LAPACK's DLASSQ accumulation keeps the convergence norm finite even when residual magnitudes differ substantially.
        double residual_scale = 0.0;
        double residual_sum_squares = 1.0;
        for (size_t active = 0; active < active_dimension; ++active) {
            const size_t coordinate = active_coordinates[active];
            const double residual = std::abs(scales[coordinate] * beta[coordinate] - average);
            if (residual == 0.0) continue;
            if (residual_scale < residual) {
                const double ratio = residual_scale / residual;
                residual_sum_squares = 1.0 + residual_sum_squares * ratio * ratio;
                residual_scale = residual;
            } else {
                const double ratio = residual / residual_scale;
                residual_sum_squares += ratio * ratio;
            }
        }
        const double standard_deviation = residual_scale * std::sqrt(residual_sum_squares / active_dimension_dbl);
        if (standard_deviation < tolerance * average) {
            converged = true;
            break;
        }

        for (size_t active_coordinate = 0; active_coordinate < active_dimension; ++active_coordinate) {
            const size_t coordinate = active_coordinates[active_coordinate];
            const double diagonal = std::abs(game(coordinate, coordinate));
            const double old_scale = scales[coordinate];
            const double c2 = (active_dimension_dbl - 1.0) * diagonal;
            const double c1 = (active_dimension_dbl - 2.0) * (beta[coordinate] - diagonal * old_scale);
            const double diagonal_scale = diagonal * old_scale;
            const double c0 = -diagonal_scale * old_scale + 2.0 * beta[coordinate] * old_scale
                              - active_dimension_dbl * average;
            const double discriminant = c1 * c1 - 4.0 * c0 * c2;
            if (!(discriminant > 0.0) || !std::isfinite(discriminant)) return safe_fallback::equilibration_invalid;

            const double denominator = c1 + std::sqrt(discriminant);
            if (denominator == 0.0 || !std::isfinite(denominator)) return safe_fallback::equilibration_invalid;
            const double new_scale = -2.0 * c0 / denominator;
            if (!(new_scale > 0.0) || !std::isfinite(new_scale)) return safe_fallback::equilibration_invalid;

            const double difference = new_scale - old_scale;
            double row_product = 0.0;
            for (size_t active_other = 0; active_other <= active_coordinate; ++active_other) {
                const size_t other = active_coordinates[active_other];
                const double magnitude = std::abs(game(coordinate, other));
                row_product += scales[other] * magnitude;
                beta[other] += difference * magnitude;
            }
            for (size_t active_other = active_coordinate + 1; active_other < active_dimension; ++active_other) {
                const size_t other = active_coordinates[active_other];
                const double magnitude = std::abs(game(other, coordinate));
                row_product += scales[other] * magnitude;
                beta[other] += difference * magnitude;
            }

            average += (row_product + beta[coordinate]) * difference / active_dimension_dbl;
            scales[coordinate] = new_scale;
        }
    }
    if (!converged) return safe_fallback::equilibration_non_convergence;

    const double common_scale = 1.0 / std::sqrt(average);
    for (size_t active = 0; active < active_dimension; ++active) {
        const size_t coordinate = active_coordinates[active];
        const double unrounded_scale = scales[coordinate] * common_scale;
        if (!(unrounded_scale > 0.0) || !std::isfinite(unrounded_scale)) return safe_fallback::equilibration_invalid;
        const double exponent = std::log2(unrounded_scale);
        if (!std::isfinite(exponent) || exponent < static_cast<double>(std::numeric_limits<int>::min())
            || exponent > static_cast<double>(std::numeric_limits<int>::max())) {
            return safe_fallback::equilibration_invalid;
        }
        scales[coordinate] = std::scalbn(1.0, static_cast<int>(exponent));
        if (!(scales[coordinate] > 0.0) || !std::isfinite(scales[coordinate])) return safe_fallback::equilibration_invalid;
    }

    for (size_t row = 0; row < dimension; ++row) {
        for (size_t column = 0; column <= row; ++column) {
            const double value = game(row, column) * (scales[row] * scales[column]);
            if (!std::isfinite(value)) return safe_fallback::equilibration_invalid;
            game(row, column) = value;
            game(column, row) = value;
        }
    }
    return safe_fallback::none;
}

void swap_active_coordinates(linalg::matrix_dbl& system, size_t active_start, size_t first, size_t second, size_t dimension,
                             size_t pivot_size)
{
    if (first == second) return;

    // Exchange one symmetric coordinate in the active lower triangle. Earlier columns contain completed factors and are left in
    // LAPACK's pivoted storage; the pivot array applies those permutations during the solve.
    for (size_t row = second + 1; row < dimension; ++row) std::swap(system(row, first), system(row, second));
    for (size_t row = first + 1; row < second; ++row) std::swap(system(row, first), system(second, row));
    std::swap(system(first, first), system(second, second));

    // A 2x2 block moves coordinate second into active_start+1 while retaining its coupling to active_start.
    if (pivot_size == 2) std::swap(system(active_start + 1, active_start), system(second, active_start));
}

/*
 * Compute P*L*D*L^T*P^T in the lower triangle of system. pivots uses LAPACK's one-based sign convention: a positive entry marks
 * a 1x1 block, while equal negative entries mark a 2x2 block. Returning false means the double result is inconclusive and exact
 * arithmetic must decide the support.
 */
bool factor_bunch_kaufman(linalg::matrix_dbl& system, size_t dimension, int* pivots)
{
    size_t k = 0;
    while (k < dimension) {
        size_t pivot_size = 1;
        size_t pivot_index = k;

        const double diagonal_magnitude = std::abs(system(k, k));
        size_t column_max_index = k;
        double column_maximum = 0.0;
        for (size_t row = k + 1; row < dimension; ++row) {
            const double magnitude = std::abs(system(row, k));
            if (magnitude > column_maximum) {
                column_maximum = magnitude;
                column_max_index = row;
            }
        }

        const double active_maximum = std::max(diagonal_magnitude, column_maximum);
        if (!std::isfinite(active_maximum) || active_maximum < kPivotCutoff) return false;

        if (diagonal_magnitude < kBunchKaufmanAlpha * column_maximum) {
            double row_maximum = 0.0;
            for (size_t column = k; column < column_max_index; ++column) {
                row_maximum = std::max(row_maximum, std::abs(system(column_max_index, column)));
            }
            for (size_t row = column_max_index + 1; row < dimension; ++row) {
                row_maximum = std::max(row_maximum, std::abs(system(row, column_max_index)));
            }
            if (!std::isfinite(row_maximum)) return false;

            if (diagonal_magnitude >= kBunchKaufmanAlpha * column_maximum * (column_maximum / row_maximum)) {
                pivot_index = k;
            } else if (std::abs(system(column_max_index, column_max_index)) >= kBunchKaufmanAlpha * row_maximum) {
                pivot_index = column_max_index;
            } else {
                pivot_index = column_max_index;
                pivot_size = 2;
            }
        }

        const size_t block_last = k + pivot_size - 1;
        swap_active_coordinates(system, k, block_last, pivot_index, dimension, pivot_size);

        if (pivot_size == 1) {
            const double pivot = system(k, k);
            if (!std::isfinite(pivot) || std::abs(pivot) < kPivotCutoff) return false;
            const double inverse_pivot = 1.0 / pivot;

            // LAPACK's rank-1 update uses the unscaled column W=L*D, then stores L in that column.
            for (size_t column = k + 1; column < dimension; ++column) {
                const double multiplier = system(column, k) * inverse_pivot;
                for (size_t row = column; row < dimension; ++row) system(row, column) -= system(row, k) * multiplier;
            }
            for (size_t row = k + 1; row < dimension; ++row) system(row, k) *= inverse_pivot;

            pivots[k] = static_cast<int>(pivot_index + 1);
        } else {
            const double off_diagonal = system(k + 1, k);
            if (!std::isfinite(off_diagonal) || std::abs(off_diagonal) < kPivotCutoff) return false;

            const double lower_diagonal = system(k + 1, k + 1) / off_diagonal;
            const double upper_diagonal = system(k, k) / off_diagonal;
            const double determinant_factor = lower_diagonal * upper_diagonal - 1.0;
            if (!std::isfinite(determinant_factor) || determinant_factor == 0.0) return false;
            const double inverse_block_scale = 1.0 / (determinant_factor * off_diagonal);

            // Each row is multiplied by inv(D), used immediately for the rank-2 Schur update, and then stored as two columns of L.
            for (size_t column = k + 2; column < dimension; ++column) {
                const double first_multiplier =
                    inverse_block_scale * (lower_diagonal * system(column, k) - system(column, k + 1));
                const double second_multiplier =
                    inverse_block_scale * (upper_diagonal * system(column, k + 1) - system(column, k));
                for (size_t row = column; row < dimension; ++row) {
                    system(row, column) -= system(row, k) * first_multiplier + system(row, k + 1) * second_multiplier;
                }
                system(column, k) = first_multiplier;
                system(column, k + 1) = second_multiplier;
            }

            pivots[k] = pivots[k + 1] = -static_cast<int>(pivot_index + 1);
        }

        k += pivot_size;
    }

    return true;
}

/* Solve the single right-hand side in place using the lower-triangle factorization above. */
bool solve_bunch_kaufman(const linalg::matrix_dbl& system, size_t dimension, const int* pivots, double* solution)
{
    // Solve L*D*y=P^T*b, applying each stored interchange immediately before its triangular update.
    size_t k = 0;
    while (k < dimension) {
        if (pivots[k] > 0) {
            const size_t pivot_index = static_cast<size_t>(pivots[k] - 1);
            if (pivot_index != k) std::swap(solution[k], solution[pivot_index]);
            for (size_t row = k + 1; row < dimension; ++row) solution[row] -= system(row, k) * solution[k];
            solution[k] /= system(k, k);
            ++k;
        } else {
            const size_t pivot_index = static_cast<size_t>(-pivots[k] - 1);
            if (pivot_index != k + 1) std::swap(solution[k + 1], solution[pivot_index]);
            for (size_t row = k + 2; row < dimension; ++row) {
                solution[row] -= system(row, k) * solution[k] + system(row, k + 1) * solution[k + 1];
            }

            const double off_diagonal = system(k + 1, k);
            const double first_diagonal = system(k, k) / off_diagonal;
            const double second_diagonal = system(k + 1, k + 1) / off_diagonal;
            const double determinant_factor = first_diagonal * second_diagonal - 1.0;
            const double first_rhs = solution[k] / off_diagonal;
            const double second_rhs = solution[k + 1] / off_diagonal;
            solution[k] = (second_diagonal * first_rhs - second_rhs) / determinant_factor;
            solution[k + 1] = (first_diagonal * second_rhs - first_rhs) / determinant_factor;
            k += 2;
        }
    }

    // Solve L^T*x=y in reverse order, then undo each corresponding interchange.
    k = dimension;
    while (k > 0) {
        const size_t block_last = k - 1;
        if (pivots[block_last] > 0) {
            for (size_t row = block_last + 1; row < dimension; ++row) {
                solution[block_last] -= system(row, block_last) * solution[row];
            }
            const size_t pivot_index = static_cast<size_t>(pivots[block_last] - 1);
            if (pivot_index != block_last) std::swap(solution[block_last], solution[pivot_index]);
            --k;
        } else {
            const size_t block_first = block_last - 1;
            for (size_t row = block_last + 1; row < dimension; ++row) {
                solution[block_last] -= system(row, block_last) * solution[row];
                solution[block_first] -= system(row, block_first) * solution[row];
            }
            const size_t pivot_index = static_cast<size_t>(-pivots[block_last] - 1);
            if (pivot_index != block_last) std::swap(solution[block_last], solution[pivot_index]);
            k -= 2;
        }
    }

    for (size_t row = 0; row < dimension; ++row) {
        if (!std::isfinite(solution[row])) return false;
    }
    return true;
}

} // namespace

find_candidate_test::find_candidate_test(const linalg::matrix_frc& game_matrix) noexcept
    : dimension_(game_matrix.rows())
{
}

void find_candidate_test::convert_game_matrix(const find_candidate_safe& safe_search)
{
    // The safe solver has already cleared one common denominator from the complete game. Normalize that integer matrix and compute
    // D*A*D once; every support reuses a principal submatrix and the matching entries of D.
    safe_fallback_ = safe_search.prepare_normalized_double_game(kPrecisionSpanCutoff, game_dbl_);
    if (safe_fallback_ == safe_fallback::none) safe_fallback_ = equilibrate_game_matrix(game_dbl_, dimension_, game_scales_.data());
}

bool find_candidate_test::find(const bitset64& support, size_t support_size)
{
    /*
     * Eliminate the normalization/payoff border exactly as find_candidate_safe does. With reference strategy m and columns
     * Z=(e_i-e_m), every normalized support vector has the form x=e_m+Z*y, and the candidate equations reduce to
     *
     *     H*y=r,       H=Z^T*A_S*Z,       r=-Z^T*A_S*e_m.
     *
     * H is symmetric and has dimension |S|-1 instead of the bordered system's |S|+1. This experimental path solves H with an
     * adapted LAPACK Bunch-Kaufman LDL^T factorization but deliberately does not use its inertia for stability.
     */
    uint8_t support_indices[bs64::kMaxBitsetDimension];
    const size_t support_count = bs64::extract_set_indices(support, dimension_, support_indices);
    bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
    assert(support_count == support_size);
    static_cast<void>(support_size);

    double solution[bs64::kMaxBitsetDimension];
    double scale_ratios[bs64::kMaxBitsetDimension];
    double reference_probability = 1.0;
    const size_t reference = support_indices[0];
    const size_t reduced_dimension = support_count - 1;
    const double inverse_reference_scale = 1.0 / game_scales_[reference];

    if (reduced_dimension > 0) {
        if (reduced_system_.rows() != reduced_dimension) reduced_system_ = linalg::matrix_dbl(reduced_dimension, reduced_dimension);

        const double scaled_reference_diagonal = game_dbl_(reference, reference);
        for (size_t position = 0; position < reduced_dimension; ++position) {
            scale_ratios[position] = game_scales_[support_indices[position + 1]] * inverse_reference_scale;
        }

        /*
         * game_dbl_ stores B=D*A*D. In the transformed coordinates x=D*z, normalization is (D*1)^T*z=1. Eliminating its
         * reference coordinate produces D_S*(Z^T*A*Z)*D_S directly from the relevant principal submatrix of B.
         */
        for (size_t row = 0; row < reduced_dimension; ++row) {
            const size_t i = support_indices[row + 1];
            const double row_ratio = scale_ratios[row];
            solution[row] = (row_ratio * scaled_reference_diagonal - game_dbl_(i, reference)) * inverse_reference_scale;
            for (size_t column = 0; column <= row; ++column) {
                const size_t j = support_indices[column + 1];
                const double column_ratio = scale_ratios[column];
                reduced_system_(row, column) = game_dbl_(i, j) - row_ratio * game_dbl_(reference, j)
                                                - column_ratio * game_dbl_(i, reference)
                                                + row_ratio * column_ratio * scaled_reference_diagonal;
            }
        }

        int pivots[bs64::kMaxBitsetDimension];
        if (!factor_bunch_kaufman(reduced_system_, reduced_dimension, pivots)
            || !solve_bunch_kaufman(reduced_system_, reduced_dimension, pivots, solution)) {
            return true;
        }

        for (size_t position = 0; position < reduced_dimension; ++position) {
            const double probability = game_scales_[support_indices[position + 1]] * solution[position];
            if (!std::isfinite(probability)) return true;
            if (probability < kProbabilityTolerance) return false;
            reference_probability -= probability;
        }
        if (!std::isfinite(reference_probability)) return true;
        if (reference_probability < kProbabilityTolerance) return false;
    }

    // B*z=D*A*x, hence division by the row's scale recovers each payoff in the globally normalized original game.
    const double reference_coordinate = reference_probability * inverse_reference_scale;
    double payoff = game_dbl_(reference, reference) * reference_coordinate;
    for (size_t position = 0; position < reduced_dimension; ++position) {
        payoff += game_dbl_(reference, support_indices[position + 1]) * solution[position];
    }
    payoff *= inverse_reference_scale;
    if (!std::isfinite(payoff)) return true;

    const double threshold = payoff + 1e-4 * static_cast<double>(dimension_);
    if (!std::isfinite(threshold)) return true;
    while (complement != 0) {
        const size_t i = bs64::find_pos_first_set_bit(complement);
        complement &= complement - 1;
        double rowsum = game_dbl_(i, reference) * reference_coordinate;
        for (size_t position = 0; position < reduced_dimension; ++position) {
            rowsum += game_dbl_(i, support_indices[position + 1]) * solution[position];
        }
        rowsum /= game_scales_[i];
        if (!std::isfinite(rowsum)) return true;
        if (rowsum > threshold) return false;
    }

    return true;
}

} // namespace candidate_search
