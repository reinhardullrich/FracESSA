#ifndef RATIONAL_LINALG_LINEAR_SOLVER_HPP
#define RATIONAL_LINALG_LINEAR_SOLVER_HPP

#include <linalg/matrix_fraction.hpp>

namespace linalg {

/*
 * Solve the exact support-equilibrium system
 *
 *     [ A_S  -1 ] [x] = [0]
 *     [ 1^T    0 ] [u]   [1].
 *
 * For support size k, M has k+1 rows and k+2 columns; its last column is the
 * right-hand side. The first k unknowns are support probabilities and the last
 * unknown u is the equilibrium payoff. A valid interior candidate requires all
 * probabilities to be strictly positive, but u may have any sign.
 *
 * M is reusable scratch storage and is intentionally destroyed by Gaussian
 * elimination. `false` means the system is singular or does not produce an
 * interior support strategy; in that case the caller must ignore x because it
 * may be incomplete. FLINT's fused submul operation avoids a temporary rational
 * in the innermost update.
 */
inline bool solve_linear_frc(matrix_frc& M, matrix_frc& x) {
    const size_t n = M.rows();

    // Exact arithmetic has no rounding-instability problem, so any nonzero pivot
    // is valid. Taking the first one avoids a full magnitude scan in every column.
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

    // Allocate the result only after elimination has shown that the system is
    // nonsingular. This avoids allocating x for singular supports rejected above.
    x = matrix_frc(n, 1);
    for (size_t i = n; i-- > 0; ) {
        fraction sum = M(i, n);

        for (size_t j = i + 1; j < n; ++j) {
            sum.submul(M(i, j), x(j, 0));
        }

        fraction::div(x(i, 0), sum, M(i, i)); 
        // The final component is u; only the preceding support weights must be positive.
        if (i < n - 1 && x(i, 0).sgn() <= 0) {
            return false; 
        }
    }

    return true;
}

} // namespace linalg

#endif // RATIONAL_LINALG_LINEAR_SOLVER_HPP
