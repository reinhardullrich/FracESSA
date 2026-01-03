#ifndef RATIONAL_LINALG_LINEAR_SOLVER_DOUBLE_HPP
#define RATIONAL_LINALG_LINEAR_SOLVER_DOUBLE_HPP

#include <rational_linalg/matrix_double.hpp>
#include <cmath>
#include <algorithm>

namespace rational_linalg {

class linear_solver_double {
public:
    explicit linear_solver_double(matrix_double& Ab) : M(Ab) {}

    bool solve(matrix_double& x) {
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
        
        x = matrix_double(n, 1);
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

private:
    matrix_double& M;
};

} // namespace rational_linalg

#endif // RATIONAL_LINALG_LINEAR_SOLVER_DOUBLE_HPP
