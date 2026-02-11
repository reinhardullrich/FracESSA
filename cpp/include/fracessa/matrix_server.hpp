#ifndef MATRIX_SERVER_HPP
#define MATRIX_SERVER_HPP

#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_double.hpp>
#include <fracessa/bitset64.hpp>
#include <cstdint>

/**
 * MatrixServer - Centralized matrix storage and operations for fracessa.
 * 
 * Responsibilities:
 *   - Manages conversion between fractional and double matrix representations.
 *   - Efficiently extracts subsystems (bordered linear systems, bee matrices).
 *   - Caches intermediate matrices to reduce heap pressure.
 */
class MatrixServer {
public:
    MatrixServer(const linalg::matrix_frc& game_matrix)
        : game_frc_(game_matrix)
    {
        dimensions_ = game_matrix.rows();
        game_dbl_ = linalg::matrix_dbl(dimensions_, dimensions_);
        auto& rational_data = game_frc_.data();
        auto& double_data = game_dbl_.data();
        for (size_t i = 0; i < rational_data.size(); ++i) {
            double_data[i] = rational_data[i].to_dbl();
        }
    }

    linalg::matrix_frc& get_linear_system_frc(const bitset64& support, size_t support_size) {
        // Build bordered linear system for support-restricted equilibrium solve:
        // [ A_S  -1 ] [x] = [0]
        // [ 1^T   0 ] [u]   [1]
        // where x are support probabilities and u is equilibrium payoff.
        // Reuse scratch matrix across supports; resize only on shape change.
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

    linalg::matrix_dbl& get_linear_system_dbl(const bitset64& support, size_t support_size) {
        // Same bordered system as above, but in double precision for fast rejection.
        // Double scratch mirror for fast candidate prefilter stage.
        if (sub_bordered_dbl_.rows() != support_size + 1) {
            sub_bordered_dbl_ = linalg::matrix_dbl(support_size + 1, support_size + 2);
        }

        uint8_t support_indices[bs64::kMaxBitsetDimension];
        bs64::extract_set_indices(support, dimensions_, support_indices);

        for (size_t ab_row = 0; ab_row < support_size; ++ab_row) {
            const size_t i = support_indices[ab_row];
            for (size_t ab_col = 0; ab_col < support_size; ++ab_col) {
                const size_t j = support_indices[ab_col];
                sub_bordered_dbl_(ab_row, ab_col) = game_dbl_(i, j);
            }
            sub_bordered_dbl_(ab_row, support_size) = -1.0;
            sub_bordered_dbl_(ab_row, support_size + 1) = 0.0;
        }
        
        for (size_t i = 0; i < support_size; ++i) {
            sub_bordered_dbl_(support_size, i) = 1.0;
        }
        sub_bordered_dbl_(support_size, support_size) = 0.0;
        sub_bordered_dbl_(support_size, support_size + 1) = 1.0;
        
        return sub_bordered_dbl_;
    }

    linalg::matrix_frc& get_bee_matrix_frc(const bitset64& extended_support_reduced, size_t m) {
        // Construct the exact Bee matrix used by ESS stability tests.
        // The formula corresponds to the quadratic form restricted to
        // perturbation directions around the candidate support.
        uint8_t reduced_indices[bs64::kMaxBitsetDimension];
        size_t size = bs64::extract_set_indices(extended_support_reduced, dimensions_, reduced_indices);
        // Bee matrix is square with size |extended_support|-1.
        if (bee_frc_.rows() != size) {
            bee_frc_ = linalg::matrix_frc(size, size);
        }

        // Precompute the constant term: 2 * game_frc_(m, m)
        fraction const_term = fraction::two() * game_frc_(m, m);

        for (size_t row = 0; row < size; ++row) {
            const size_t i = reduced_indices[row];
            for (size_t column = 0; column <= row; ++column) {
                const size_t j = reduced_indices[column];
                fraction& val = bee_frc_(row, column);
                // val = game_frc_(m, j) + game_frc_(j, m) + game_frc_(i, m) + game_frc_(m, i) - game_frc_(i, j) - game_frc_(j, i) - const_term
                fraction::add(val, game_frc_(m, j), game_frc_(j, m));
                val.add_inplace(game_frc_(i, m));
                val.add_inplace(game_frc_(m, i));
                val.sub_inplace(game_frc_(i, j));
                val.sub_inplace(game_frc_(j, i));
                val.sub_inplace(const_term);

                if (row != column) {
                    bee_frc_(column, row) = val;
                }
            }
        }
        return bee_frc_;
    }

    linalg::matrix_dbl& get_bee_matrix_dbl(const bitset64& extended_support_reduced, size_t m) {
        // Double mirror of Bee construction: used as a quick positive-definite filter
        // before paying the cost of exact arithmetic.
        uint8_t reduced_indices[bs64::kMaxBitsetDimension];
        size_t size = bs64::extract_set_indices(extended_support_reduced, dimensions_, reduced_indices);
        if (bee_dbl_.rows() != size) {
            bee_dbl_ = linalg::matrix_dbl(size, size);
        }

        for (size_t row = 0; row < size; ++row) {
            const size_t i = reduced_indices[row];
            for (size_t column = 0; column <= row; ++column) {
                const size_t j = reduced_indices[column];
                bee_dbl_(row, column) = bee_dbl_(column, row) =
                    game_dbl_(m, j) + game_dbl_(j, m) + game_dbl_(i, m) + game_dbl_(m, i) -
                    game_dbl_(i, j) - game_dbl_(j, i) - 2.0 * game_dbl_(m, m);
            }
        }
        return bee_dbl_;
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
    linalg::matrix_dbl bee_dbl_;
    size_t dimensions_;

};

#endif // MATRIX_SERVER_HPP
