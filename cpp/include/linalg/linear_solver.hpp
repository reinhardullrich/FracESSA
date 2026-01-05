#ifndef RATIONAL_LINALG_LINEAR_SOLVER_HPP
#define RATIONAL_LINALG_LINEAR_SOLVER_HPP

#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_double.hpp>
#include <cmath>
#include <algorithm>
#include <vector>

namespace linalg {

/**
 * solve_linear_dbl - Optimized Standard Gaussian elimination for double matrices.
 * 
 * Side Effects:
 *   - Modifies the input matrix M IN-PLACE (Forward Elimination).
 *   - M will be in an upper triangular form after this call.
 * 
 * Performance:
 *   - Uses physical row swaps for better cache locality.
 */
inline bool solve_linear_dbl(matrix_dbl& M, matrix_dbl& x) {
    const size_t n = M.rows();
    
    for (size_t k = 0; k < n - 1; ++k) {
        size_t max_row = k;
        double max_val = std::abs(M(k, k));
        
        for (size_t i = k + 1; i < n; ++i) {
            double val = std::abs(M(i, k));
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        
        if (max_val < 1e-12) return false;
        
        if (max_row != k) {
            M.swap_rows(k, max_row);
        }
        
        const double pivot = M(k, k);
        for (size_t i = k + 1; i < n; ++i) {
            const double factor = M(i, k) / pivot;
            for (size_t j = k + 1; j <= n; ++j) {
                M(i, j) -= factor * M(k, j);
            }
            M(i, k) = 0.0;
        }
    }
    
    if (std::abs(M(n - 1, n - 1)) < 1e-12) return false;
    
    x = matrix_dbl(n, 1);
    for (size_t i = n; i-- > 0; ) {
        double sum = M(i, n);
        for (size_t j = i + 1; j < n; ++j) {
            sum -= M(i, j) * x(j, 0);
        }
        
        const double pivot = M(i, i);
        double temp_x = sum / pivot;
        
        if (temp_x < -1e-10) return false; // Early exit if not strictly positive
        
        x(i, 0) = temp_x;
    }
    
    return true;
}

/**
 * solve_linear_frc - Optimized Gaussian elimination for rational matrices.
 * 
 * Side Effects:
 *   - Modifies the input matrix M IN-PLACE (Forward Elimination).
 * 
 * Optimizations:
 *   - Uses fmpq_submul (FLINT) for direct subtraction-multiplication without temporaries.
 *   - Uses physical row swaps instead of virtual pivoting for better cache locality.
 */
inline bool solve_linear_frc(matrix_frc& M, matrix_frc& x) {
    const size_t n = M.rows();

    // Forward elimination
    for (size_t k = 0; k < n - 1; ++k) {
        size_t max_row = k;
        
        if (M(k, k).is_zero()) {
            bool found = false;
            for (size_t i = k + 1; i < n; ++i) {
                if (!M(i, k).is_zero()) {
                    max_row = i;
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        
        if (max_row != k) {
            M.swap_rows(k, max_row);
        }
        
        const fraction& pivot = M(k, k);

        for (size_t i = k + 1; i < n; ++i) {
            if (M(i, k).is_zero()) continue;

            fraction factor;
            fraction::div(factor, M(i, k), pivot);

            for (size_t j = k + 1; j <= n; ++j) {
                M(i, j).submul(factor, M(k, j));
            }
            
            M(i, k) = fraction::zero();
        }
    }

    if (M(n - 1, n - 1).is_zero()) {
        return false;
    }

    // Back substitution
    x = matrix_frc(n, 1);
    for (size_t i = n; i-- > 0; ) {
        fraction sum = M(i, n);

        for (size_t j = i + 1; j < n; ++j) {
            sum.submul(M(i, j), x(j, 0));
        }

        fraction::div(x(i, 0), sum, M(i, i)); // Early exit: x must be > 0.
        // Full support for this subgame <=> x is in the interior of the (sub)simplex 
        // <=> all strategy probabilities x_i are strictly positive.
        if (x(i, 0).sgn() <= 0) {
            return false; 
        }
    }

    return true;
}

} // namespace linalg

#endif // RATIONAL_LINALG_LINEAR_SOLVER_HPP
