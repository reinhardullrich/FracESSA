#ifndef MATRIX_SERVER_HPP
#define MATRIX_SERVER_HPP

#include <rational_linalg/matrix_fraction.hpp>
#include <rational_linalg/matrix_double.hpp>
#include <rational_linalg/positive_definite_fraction.hpp>
#include <rational_linalg/positive_definite_double.hpp>
#include <fracessa/bitset64.hpp>

/**
 * MatrixServer - Centralized matrix storage and operations for fracessa
 */
class MatrixServer {
public:
    MatrixServer(const rational_linalg::matrix_fraction& game_matrix)
        : game_rational_(game_matrix)
    {
        dimensions_ = game_matrix.rows();
        game_double_ = rational_linalg::matrix_double(dimensions_, dimensions_);
        auto& rational_data = game_rational_.data();
        auto& double_data = game_double_.data();
        for (size_t i = 0; i < rational_data.size(); ++i) {
            double_data[i] = rational_data[i].to_double();
        }
    }

    rational_linalg::matrix_fraction& get_linear_system_fraction(const bitset64& support, size_t support_size) {
        if (subgame_augmented_rational_.rows() != support_size + 1) {
            subgame_augmented_rational_ = rational_linalg::matrix_fraction(support_size + 1, support_size + 2);
        }

        size_t ab_row = 0;
        for (size_t i = 0; i < dimensions_; ++i) {
            if (bs64::is_set_at_pos(support, i)) {
                size_t ab_col = 0;
                for (size_t j = 0; j < dimensions_; ++j) {
                    if (bs64::is_set_at_pos(support, j)) {
                        subgame_augmented_rational_(ab_row, ab_col) = game_rational_(i, j);
                        ab_col++;
                    }
                }
                subgame_augmented_rational_(ab_row, support_size) = fraction::neg_one();
                subgame_augmented_rational_(ab_row, support_size + 1) = fraction::zero();
                ab_row++;
            }
        }
        
        for (size_t i = 0; i < support_size; ++i) {
            subgame_augmented_rational_(support_size, i) = fraction::one();
        }
        subgame_augmented_rational_(support_size, support_size) = fraction::zero();
        subgame_augmented_rational_(support_size, support_size + 1) = fraction::one();
        
        return subgame_augmented_rational_;
    }

    rational_linalg::matrix_double& get_linear_system_double(const bitset64& support, size_t support_size) {
        if (subgame_augmented_double_.rows() != support_size + 1) {
            subgame_augmented_double_ = rational_linalg::matrix_double(support_size + 1, support_size + 2);
        }

        size_t ab_row = 0;
        for (size_t i = 0; i < dimensions_; ++i) {
            if (bs64::is_set_at_pos(support, i)) {
                size_t ab_col = 0;
                for (size_t j = 0; j < dimensions_; ++j) {
                    if (bs64::is_set_at_pos(support, j)) {
                        subgame_augmented_double_(ab_row, ab_col) = game_double_(i, j);
                        ab_col++;
                    }
                }
                subgame_augmented_double_(ab_row, support_size) = -1.0;
                subgame_augmented_double_(ab_row, support_size + 1) = 0.0;
                ab_row++;
            }
        }
        
        for (size_t i = 0; i < support_size; ++i) {
            subgame_augmented_double_(support_size, i) = 1.0;
        }
        subgame_augmented_double_(support_size, support_size) = 0.0;
        subgame_augmented_double_(support_size, support_size + 1) = 1.0;
        
        return subgame_augmented_double_;
    }

    rational_linalg::matrix_fraction& get_bee_matrix_fraction(const bitset64& extended_support_reduced, size_t m) {
        size_t size = bs64::count_set_bits(extended_support_reduced);
        if (bee_rational_.rows() != size) {
            bee_rational_ = rational_linalg::matrix_fraction(size, size);
        }

        size_t row = 0;
        for (size_t i = 0; i < dimensions_; i++) {
            if (bs64::is_set_at_pos(extended_support_reduced, i)) {
                size_t column = 0;
                for (size_t j = 0; j < i + 1; j++) {
                    if (bs64::is_set_at_pos(extended_support_reduced, j)) {
                        bee_rational_(row, column) = bee_rational_(column, row) = 
                            game_rational_(m, j) + game_rational_(j, m) + game_rational_(i, m) + game_rational_(m, i) -
                            game_rational_(i, j) - game_rational_(j, i) - fraction::two() * game_rational_(m, m);
                        column += 1;
                    }
                }
                row += 1;
            }
        }
        return bee_rational_;
    }

    rational_linalg::matrix_double& get_bee_matrix_double(const bitset64& extended_support_reduced, size_t m) {
        size_t size = bs64::count_set_bits(extended_support_reduced);
        if (bee_double_.rows() != size) {
            bee_double_ = rational_linalg::matrix_double(size, size);
        }

        size_t row = 0;
        for (size_t i = 0; i < dimensions_; i++) {
            if (bs64::is_set_at_pos(extended_support_reduced, i)) {
                size_t column = 0;
                for (size_t j = 0; j < i + 1; j++) {
                    if (bs64::is_set_at_pos(extended_support_reduced, j)) {
                        bee_double_(row, column) = bee_double_(column, row) = 
                            game_double_(m, j) + game_double_(j, m) + game_double_(i, m) + game_double_(m, i) -
                            game_double_(i, j) - game_double_(j, i) - 2.0 * game_double_(m, m);
                        column += 1;
                    }
                }
                row += 1;
            }
        }
        return bee_double_;
    }

    const rational_linalg::matrix_fraction& get_game_matrix_fraction() const noexcept {
        return game_rational_;
    }

    const rational_linalg::matrix_double& get_game_matrix_double() const noexcept {
        return game_double_;
    }

private:
    rational_linalg::matrix_fraction game_rational_;
    rational_linalg::matrix_double game_double_;
    rational_linalg::matrix_fraction subgame_augmented_rational_;
    rational_linalg::matrix_double subgame_augmented_double_;
    rational_linalg::matrix_fraction bee_rational_;
    rational_linalg::matrix_double bee_double_;
    size_t dimensions_;
};

#endif // MATRIX_SERVER_HPP
