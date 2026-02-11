#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/linear_solver.hpp>
#include <cstdlib>
#include <cstdint>

namespace {

inline void extract_support_partition(
    bitset64 support,
    size_t dimension,
    uint8_t (&support_indices)[bs64::kMaxBitsetDimension],
    size_t& support_count,
    uint8_t (&non_support_indices)[bs64::kMaxBitsetDimension],
    size_t& non_support_count) noexcept
{
    support_count = bs64::extract_set_indices(support, dimension, support_indices);
    non_support_count = 0;

    size_t support_pos = 0;
    for (size_t i = 0; i < dimension; ++i) {
        if (support_pos < support_count && i == static_cast<size_t>(support_indices[support_pos])) {
            ++support_pos;
        } else {
            non_support_indices[non_support_count++] = static_cast<uint8_t>(i);
        }
    }
}

} // namespace

bool fracessa::find_candidate_dbl(const bitset64& support, size_t support_size)
{
    const auto& game_matrix = matrix_server_.get_game_matrix_dbl();
    auto& Ab = matrix_server_.get_linear_system_dbl(support, support_size);
    
    linalg::matrix_dbl solution;
    bool solved = linalg::solve_linear_dbl(Ab, solution);
    
    if (!solved) return false;

    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    size_t support_count = 0;
    size_t non_support_count = 0;
    extract_support_partition(support, dimension_, support_indices, support_count, non_support_indices, non_support_count);
    double full_solution[bs64::kMaxBitsetDimension];
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        full_solution[non_support_indices[i_pos]] = 0.0;
    }
    for (size_t i_pos = 0; i_pos < support_count; ++i_pos) {
        full_solution[support_indices[i_pos]] = solution(i_pos, 0);
    }

    // --- FLOATING POINT FILTER PATH ---
    // Use a large margin (threshold) to filter out clearly invalid supports early
    // without the cost of exact arithmetic. The 1e-4 * dimension_ factor provides
    // enough buffer to avoid false negatives due to precision loss.
    const double payoff = solution(support_size, 0);
    const double threshold = payoff + 1e-4 * dimension_;
    
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        double rowsum = 0.0;
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            const size_t j = support_indices[j_pos];
            rowsum += game_matrix(i, j) * full_solution[j];
        }
        if (rowsum > threshold) return false;
    }
    return true;
}

bool fracessa::find_candidate_frc(const bitset64& support, size_t support_size)
{
    const auto& game_matrix = matrix_server_.get_game_matrix_frc();
    auto& Ab = matrix_server_.get_linear_system_frc(support, support_size);
    
    linalg::matrix_frc solution;
    bool solved = linalg::solve_linear_frc(Ab, solution);
    
    if (!solved) return false;

    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    size_t support_count = 0;
    size_t non_support_count = 0;
    extract_support_partition(support, dimension_, support_indices, support_count, non_support_indices, non_support_count);

    const fraction payoff = solution(support_size, 0);
    candidate_.extended_support = support;

    if (candidate_.vector.rows() != dimension_ || candidate_.vector.cols() != 1) {
        candidate_.vector = linalg::matrix_frc(dimension_, 1);
    }
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        candidate_.vector(non_support_indices[i_pos], 0) = fraction::zero();
    }
    for (size_t i_pos = 0; i_pos < support_count; ++i_pos) {
        candidate_.vector(support_indices[i_pos], 0) = solution(i_pos, 0);
    }
    
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        fraction rowsum = fraction::zero();
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            const size_t j = support_indices[j_pos];
            rowsum.addmul(game_matrix(i, j), solution(j_pos, 0));
        }
        
        if (rowsum > payoff) return false;
        
        if (rowsum == payoff)
            candidate_.extended_support = bs64::set_bit_at_pos(candidate_.extended_support, i);
    }

    candidate_.payoff = payoff;
    candidate_.payoff_dbl = payoff.to_dbl();        
    candidate_.extended_support_size = bs64::count_set_bits(candidate_.extended_support);
    
    return true;
}
