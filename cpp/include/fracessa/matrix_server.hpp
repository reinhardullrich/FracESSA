#ifndef MATRIX_SERVER_HPP
#define MATRIX_SERVER_HPP

#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_double.hpp>
#include <fracessa/bitset64.hpp>
#include <cmath>
#include <cstdint>

/*
 * Owns the payoff matrix and the small work matrices derived from it.
 *
 * The same shapes occur for many supports. Keeping the bordered system and Bee
 * matrix here lets the search reuse their allocations instead of asking the
 * heap for another small matrix millions of times. The exact game remains the
 * source of truth; the double copy exists only for the unsafe rejection filter.
 */
class MatrixServer {
public:
    MatrixServer(const linalg::matrix_frc& game_matrix)
        : game_frc_(game_matrix), dimensions_(game_matrix.rows()) {}

    /*
     * Build one normalized binary64 copy for the unsafe filter:
     *
     *     A' = (A - c) / s,
     *     c = one exact matrix entry,
     *     s = max |A_ij - c|.
     *
     * Adding a common constant and multiplying by a positive constant preserve
     * all best-response comparisons and ESS decisions. Subtraction and division
     * happen exactly before conversion, removing irrelevant common offset and
     * scale before binary64 rounding is introduced.
     *
     * False means the double copy is unusable (for example, the game is
     * constant or a nonzero value underflows). The analyzer then uses exact
     * candidate solving for every support.
     */
    bool initialize_unsafe_filter() {
        const auto& rational_data = game_frc_.data();
        const fraction& translation = rational_data.front();
        fraction scale;
        fraction difference;

        for (const auto& value : rational_data) {
            fraction::sub(difference, value, translation);
            if (difference.sgn() < 0) {
                fmpq_neg(difference.data(), difference.data());
            }
            if (difference > scale) {
                scale = difference;
            }
        }
        if (scale.is_zero()) {
            return false;
        }

        game_dbl_ = linalg::matrix_dbl(dimensions_, dimensions_);
        auto& double_data = game_dbl_.data();
        fraction normalized;
        for (size_t i = 0; i < rational_data.size(); ++i) {
            fraction::sub(difference, rational_data[i], translation);
            fraction::div(normalized, difference, scale);
            const double converted = normalized.to_dbl();
            if (!std::isfinite(converted) ||
                (converted == 0.0 && !normalized.is_zero()) ||
                std::fpclassify(converted) == FP_SUBNORMAL) {
                return false;
            }
            double_data[i] = converted;
        }
        return true;
    }

    linalg::matrix_frc& get_linear_system_frc(const bitset64& support, size_t support_size) {
        /*
         * For a support S of size k, build the augmented system
         *
         *     [ A_S  -1 | 0 ]
         *     [ 1^T    0 | 1 ].
         *
         * Its k+1 unknowns are the k positive support probabilities x and the
         * common payoff u. The last column is the right-hand side. Reuse the
         * (k+1)-by-(k+2) storage whenever the next support has the same size.
         */
        if (sub_bordered_frc_.rows() != support_size + 1) {
            sub_bordered_frc_ = linalg::matrix_frc(support_size + 1, support_size + 2);
        }

        uint8_t support_indices[bs64::kMaxBitsetDimension];
        bs64::extract_set_indices(support, dimensions_, support_indices);

        for (size_t ab_row = 0; ab_row < support_size; ++ab_row) {
            const size_t i = support_indices[ab_row];
            for (size_t ab_col = 0; ab_col < support_size; ++ab_col) {
                const size_t j = support_indices[ab_col];
                sub_bordered_frc_(ab_row, ab_col) = game_frc_(i, j);
            }
            sub_bordered_frc_(ab_row, support_size) = fraction::neg_one();
            sub_bordered_frc_(ab_row, support_size + 1) = fraction::zero();
        }
        
        for (size_t i = 0; i < support_size; ++i) {
            sub_bordered_frc_(support_size, i) = fraction::one();
        }
        sub_bordered_frc_(support_size, support_size) = fraction::zero();
        sub_bordered_frc_(support_size, support_size + 1) = fraction::one();
        
        return sub_bordered_frc_;
    }

    linalg::matrix_dbl& get_linear_system_dbl(size_t support_size) {
        // The unsafe solver builds the same augmented system directly from the
        // normalized double game, so this method only supplies reusable storage.
        if (sub_bordered_dbl_.rows() != support_size + 1) {
            sub_bordered_dbl_ = linalg::matrix_dbl(support_size + 1, support_size + 2);
        }
        return sub_bordered_dbl_;
    }

    linalg::matrix_frc& get_bee_matrix_frc(const bitset64& extended_support_reduced, size_t m) {
        /*
         * Construct Bomze's exact stability matrix B (`Bee` in code identifiers),
         * defined in Bomze (1992), Theorem 3.2.
         *
         * Let J be the candidate's extended support, choose m in its support,
         * and index this matrix by J without m. For i,j in that reduced set,
         *
         * B_ij = a_mj + a_jm + a_im + a_mi
         *        - a_ij - a_ji - 2*a_mm.
         *
         * The ESS condition becomes positivity of y^T B y on a cone determined
         * by the unused best replies. `check_stability()` performs that test.
         */
        uint8_t reduced_indices[bs64::kMaxBitsetDimension];
        size_t size = bs64::extract_set_indices(extended_support_reduced, dimensions_, reduced_indices);
        // B is square with one row and column for each member of J except m.
        if (bee_frc_.rows() != size) {
            bee_frc_ = linalg::matrix_frc(size, size);
        }

        // The final term is shared by every entry.
        fraction const_term = fraction::two() * game_frc_(m, m);

        for (size_t row = 0; row < size; ++row) {
            const size_t i = reduced_indices[row];
            for (size_t column = 0; column <= row; ++column) {
                const size_t j = reduced_indices[column];
                fraction& val = bee_frc_(row, column);
                // Compute only one triangle, then mirror it because B is symmetric.
                fraction::add(val, game_frc_(m, j), game_frc_(j, m));
                val += game_frc_(i, m);
                val += game_frc_(m, i);
                val -= game_frc_(i, j);
                val -= game_frc_(j, i);
                val -= const_term;

                if (row != column) {
                    bee_frc_(column, row) = val;
                }
            }
        }
        return bee_frc_;
    }

    const linalg::matrix_frc& get_game_matrix_frc() const noexcept {
        return game_frc_;
    }

    const linalg::matrix_dbl& get_game_matrix_dbl() const noexcept {
        return game_dbl_;
    }

private:
    linalg::matrix_frc game_frc_;
    linalg::matrix_dbl game_dbl_;
    linalg::matrix_frc sub_bordered_frc_;
    linalg::matrix_dbl sub_bordered_dbl_;
    linalg::matrix_frc bee_frc_;
    size_t dimensions_;

};

#endif // MATRIX_SERVER_HPP
