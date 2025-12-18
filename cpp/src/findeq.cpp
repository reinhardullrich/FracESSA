#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/linear_solver.hpp>
#include <cstdlib>
#include <type_traits>

// Templated function for all types (double, fraction)
template<typename T>
bool fracessa::find_candidate(const bitset64& support, size_t support_size)
{
    const auto& game_matrix = matrix_server_.get_game_matrix<T>();

    // Use LinearSolver for all types
    // get_linear_system returns a reference to the matrix in MatrixServer
    auto& Ab = matrix_server_.get_linear_system<T>(support, support_size);
    
    // The optimized LinearSolver (GaussRational/GaussDouble) takes a reference
    // and modifies the matrix in-place for efficiency.
    rational_linalg::LinearSolver<T> solver(Ab);
    rational_linalg::Matrix<T> solution;
    bool solved = solver.solve(solution);
    
    if (!solved)
        return false;    

    // Build full solution vector from support by padding with zeros
    rational_linalg::Matrix<T> solution_full_n = rational_linalg::Matrix<T>(dimension_, 1);    
    size_t tracker = 0;
    for (size_t i = 0; i < dimension_; i++) {
        if (bs64::is_set_at_pos(support, i)) {
            solution_full_n(i, 0) = solution(tracker, 0);
            tracker += 1;
        } else 
            solution_full_n(i, 0) = T(0);
    }

    // Check (Ap)_i <= v for all rows i not in the support
    const T payoff = solution(support_size, 0);
    candidate_.extended_support = support; // reset/init extended support
    
    // Split logic completely for double vs fraction
    if constexpr (std::is_same_v<T, double>) {
        // --- FLOATING POINT FILTER PATH ---
        const double threshold = payoff + 1e-4 * dimension_; // huge margin to avoid false negatives
        
        for (size_t i = 0; i < dimension_; i++) {
            if (!bs64::is_set_at_pos(support, i)) { // Row not in the support
                double rowsum = 0.0;
                for (size_t j = 0; j < dimension_; j++) {
                    if (bs64::is_set_at_pos(support, j)) { // Column is in the support
                        rowsum += game_matrix(i,j) * solution_full_n(j, 0);
                    }
                }
                
                if (rowsum > threshold) 
                    return false;
            }
        }
    } else {
        // --- EXACT ARITHMETIC PATH ---
        for (size_t i = 0; i < dimension_; i++) {
            if (!bs64::is_set_at_pos(support, i)) { // Row not in the support
                T rowsum = T(0);
                for (size_t j = 0; j < dimension_; j++) {
                    if (bs64::is_set_at_pos(support, j)) { // Column is in the support
                        // Optimized FLINT arithmetic: rowsum += matrix(i,j) * solution(j)
                        fmpq_addmul(rowsum.data(), game_matrix(i,j).data(), solution_full_n(j, 0).data());
                    }
                }
                
                if (rowsum > payoff)
                    return false;
                
                // Only track extended support in the exact phase
                if (rowsum == payoff)
                    candidate_.extended_support = bs64::set_bit_at_pos(candidate_.extended_support, i);
            }
        }

        // Store candidate data (Only needed when we have the final exact candidate)
        candidate_.vector = solution_full_n;
        candidate_.payoff = payoff;
        candidate_.payoff_double = payoff.to_double();        
        candidate_.extended_support_size = bs64::count_set_bits(candidate_.extended_support);
    }
    
    return true;
}

// Explicit template instantiations
template bool fracessa::find_candidate<double>(const bitset64&, size_t);
template bool fracessa::find_candidate<fraction>(const bitset64&, size_t);