#include <fracessa/find_candidate_exact.hpp>

#include <cstdint>

#include <linalg/linear_solver.hpp>

namespace candidate_search {

/*
 * Exact candidate construction for a proposed support S.
 *
 * A mixed strategy x with support S is a symmetric Nash equilibrium when
 *
 *   A_S x = u*1,          every used strategy earns the same payoff u,
 *   1^T x = 1,            probabilities sum to one,
 *   x_i > 0 for i in S,   S is the actual support,
 *   (A x)_i <= u outside S.
 */
find_candidate_exact::find_candidate_exact(const linalg::matrix_frc& game_matrix) noexcept
    : game_frc_(game_matrix)
    , dimension_(game_matrix.rows())
{
}

bool find_candidate_exact::find(
    const bitset64& support,
    size_t support_size,
    candidate& result)
{
    const size_t system_dimension = support_size + 1;
    if (linear_system_.rows() != system_dimension) {
        linear_system_ =
            linalg::matrix_frc(system_dimension, system_dimension + 1);
    }

    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    const size_t support_count =
        bs64::extract_set_indices(support, dimension_, support_indices);
    const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
    const size_t non_support_count =
        bs64::extract_set_indices(complement, dimension_, non_support_indices);

    for (size_t row = 0; row < support_size; ++row) {
        const size_t i = support_indices[row];
        for (size_t column = 0; column < support_size; ++column) {
            linear_system_(row, column) =
                game_frc_(i, support_indices[column]);
        }
        linear_system_(row, support_size) = fraction::neg_one();
        linear_system_(row, system_dimension) = fraction::zero();
    }
    for (size_t column = 0; column < support_size; ++column) {
        linear_system_(support_size, column) = fraction::one();
    }
    linear_system_(support_size, support_size) = fraction::zero();
    linear_system_(support_size, system_dimension) = fraction::one();

    linalg::matrix_frc solution;
    if (!linalg::solve_linear_frc(linear_system_, solution)) return false;

    const fraction payoff = solution(support_size, 0);
    result.extended_support = support;

    if (result.vector.rows() != dimension_ || result.vector.cols() != 1) {
        result.vector = linalg::matrix_frc(dimension_, 1);
    }
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        result.vector(non_support_indices[i_pos], 0) = fraction::zero();
    }
    for (size_t i_pos = 0; i_pos < support_count; ++i_pos) {
        result.vector(support_indices[i_pos], 0) = solution(i_pos, 0);
    }

    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        fraction rowsum = fraction::zero();
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            rowsum.addmul(
                game_frc_(i, support_indices[j_pos]), solution(j_pos, 0));
        }

        if (rowsum > payoff) return false;
        if (rowsum == payoff) {
            result.extended_support =
                bs64::set_bit_at_pos(result.extended_support, i);
        }
    }

    result.payoff = payoff;
    result.payoff_dbl = payoff.to_dbl();
    result.extended_support_size =
        bs64::count_set_bits(result.extended_support);
    return true;
}

} // namespace candidate_search
