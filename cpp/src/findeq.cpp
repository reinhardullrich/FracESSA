#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/linear_solver_fraction.hpp>
#include <rational_linalg/linear_solver_double.hpp>
#include <cstdlib>

bool fracessa::find_candidate_double(const bitset64& support, size_t support_size)
{
    const auto& game_matrix = matrix_server_.get_game_matrix_double();
    auto& Ab = matrix_server_.get_linear_system_double(support, support_size);
    
    rational_linalg::linear_solver_double solver(Ab);
    rational_linalg::matrix_double solution;
    bool solved = solver.solve(solution);
    
    if (!solved) return false;

    rational_linalg::matrix_double solution_full_n(dimension_, 1);
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (bs64::is_set_at_pos(support, i)) {
            solution_full_n(i, 0) = solution(tracker, 0);
            tracker += 1;
        } else {
            solution_full_n(i, 0) = 0.0;
        }
    }

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

bool fracessa::find_candidate_fraction(const bitset64& support, size_t support_size)
{
    const auto& game_matrix = matrix_server_.get_game_matrix_fraction();
    auto& Ab = matrix_server_.get_linear_system_fraction(support, support_size);
    
    rational_linalg::linear_solver_fraction solver(Ab);
    rational_linalg::matrix_fraction solution;
    bool solved = solver.solve(solution);
    
    if (!solved) return false;

    rational_linalg::matrix_fraction solution_full_n(dimension_, 1);
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
    candidate_.payoff_double = payoff.to_double();        
    candidate_.extended_support_size = bs64::count_set_bits(candidate_.extended_support);
    
    return true;
}