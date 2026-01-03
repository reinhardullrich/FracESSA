#ifndef RATIONAL_LINALG_LINEAR_SOLVER_FRACTION_HPP
#define RATIONAL_LINALG_LINEAR_SOLVER_FRACTION_HPP

#include <rational_linalg/matrix_fraction.hpp>
#include <vector>

namespace rational_linalg {

class linear_solver_fraction {
public:
    explicit linear_solver_fraction(matrix_fraction& Ab) : M(Ab) {}

    bool solve(matrix_fraction& x) {
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
                
                M(i, k).set_zero();
            }
        }

        if (M(n - 1, n - 1).is_zero()) {
            return false;
        }

        // Back substitution
        x = matrix_fraction(n, 1);
        for (size_t i = n; i-- > 0; ) {
            fraction sum = M(i, n);

            for (size_t j = i + 1; j < n; ++j) {
                sum.submul(M(i, j), x(j, 0));
            }

            fraction::div(x(i, 0), sum, M(i, i));
            
            if (x(i, 0).sgn() <= 0) {
                return false; 
            }
        }

        return true;
    }

private:
    matrix_fraction& M;
};

} // namespace rational_linalg

#endif // RATIONAL_LINALG_LINEAR_SOLVER_FRACTION_HPP
