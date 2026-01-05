#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/linear_solver.hpp>
#include <cstdlib>

bool fracessa::find_candidate_dbl(const bitset64& support, size_t support_size)
{
    const auto& game_matrix = matrix_server_.get_game_matrix_dbl();
    auto& Ab = matrix_server_.get_linear_system_dbl(support, support_size);
    
    linalg::matrix_dbl solution;
    bool solved = linalg::solve_linear_dbl(Ab, solution);
    
    if (!solved) return false;

    linalg::matrix_dbl solution_full_n(dimension_, 1);
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (bs64::is_set_at_pos(support, i)) {
            solution_full_n(i, 0) = solution(tracker, 0);
            tracker += 1;
        } else {
            solution_full_n(i, 0) = 0.0;
        }
    }

    // --- FLOATING POINT FILTER PATH ---
    // Use a large margin (threshold) to filter out clearly invalid supports early
    // without the cost of exact arithmetic. The 1e-4 * dimension_ factor provides
    // enough buffer to avoid false negatives due to precision loss.
    const double payoff = solution(support_size, 0);
    const double threshold = payoff + 1e-4 * dimension_;
    
    for (size_t i = 0; i < dimension_; i++) {
        if (!bs64::is_set_at_pos(support, i)) {
            double rowsum = 0.0;
            for (size_t j = 0; j < dimension_; j++) {
                if (bs64::is_set_at_pos(support, j)) {
                    rowsum += game_matrix(i, j) * solution_full_n(j, 0);
                }
            }
            if (rowsum > threshold) return false;
        }
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

    linalg::matrix_frc solution_full_n(dimension_, 1);
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (bs64::is_set_at_pos(support, i)) {
            solution_full_n(i, 0) = solution(tracker, 0);
            tracker += 1;
        } else {
            solution_full_n(i, 0) = fraction::zero();
        }
    }

    const fraction payoff = solution(support_size, 0);
    candidate_.extended_support = support;
    
    for (size_t i = 0; i < dimension_; i++) {
        if (!bs64::is_set_at_pos(support, i)) {
            fraction rowsum = fraction::zero();
            for (size_t j = 0; j < dimension_; j++) {
                if (bs64::is_set_at_pos(support, j)) {
                    rowsum.addmul(game_matrix(i, j), solution_full_n(j, 0));
                }
            }
            
            if (rowsum > payoff) return false;
            
            if (rowsum == payoff)
                candidate_.extended_support = bs64::set_bit_at_pos(candidate_.extended_support, i);
        }
    }

    candidate_.vector = solution_full_n;
    candidate_.payoff = payoff;
    candidate_.payoff_dbl = payoff.to_dbl();        
    candidate_.extended_support_size = bs64::count_set_bits(candidate_.extended_support);
    
    return true;
}