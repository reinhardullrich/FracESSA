#ifndef RATIONAL_LINALG_POSITIVE_DEFINITE_DOUBLE_HPP
#define RATIONAL_LINALG_POSITIVE_DEFINITE_DOUBLE_HPP

#include <rational_linalg/matrix_double.hpp>
#include <limits>
#include <cmath>

namespace rational_linalg {

inline bool matrix_double::is_positive_definite() const {
    const size_t n = rows_;
    if (n == 0) return true;
    if (n != cols_) return false;
    
    matrix_double L(n, n);
    const double tolerance = std::numeric_limits<double>::epsilon() * infinity_norm();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < j; ++k) {
                sum += L(i, k) * L(j, k);
            }
            L(i, j) = (1.0 / L(j, j)) * ((*this)(i, j) - sum);
        }
        
        double sum = 0.0;
        for (size_t k = 0; k < i; ++k) {
            sum += L(i, k) * L(i, k);
        }
        double diagonal_value = (*this)(i, i) - sum;
        
        if (diagonal_value <= tolerance) {
            return false;
        }
        
        L(i, i) = std::sqrt(diagonal_value);
    }
    
    return true;
}

} // namespace rational_linalg

#endif // RATIONAL_LINALG_POSITIVE_DEFINITE_DOUBLE_HPP
