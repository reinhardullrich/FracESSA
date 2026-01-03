#ifndef RATIONAL_LINALG_ADJUGATE_FRACTION_HPP
#define RATIONAL_LINALG_ADJUGATE_FRACTION_HPP

#include <rational_linalg/matrix_fraction.hpp>
#include <rational_linalg/lu_factor_fraction.hpp>

namespace rational_linalg {

inline matrix_fraction adjugate(const matrix_fraction& A) {
    const size_t n = A.rows();
    if (n == 0) return matrix_fraction(0, 0);
    if (n == 1) {
        matrix_fraction result(1, 1);
        result(0, 0) = fraction::one();
        return result;
    }
    
    matrix_fraction cofactor_matrix(n, n);
    matrix_fraction minor(n - 1, n - 1);
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            size_t minor_row = 0;
            for (size_t row = 0; row < n; ++row) {
                if (row == i) continue;
                size_t minor_col = 0;
                for (size_t col = 0; col < n; ++col) {
                    if (col == j) continue;
                    minor(minor_row, minor_col) = A(row, col);
                    ++minor_col;
                }
                ++minor_row;
            }
            
            lu_factor_fraction lu(minor);
            fraction det_minor = lu.determinant();
            fraction sign = ((i + j) & 1) == 0 ? fraction::one() : fraction::neg_one();
            cofactor_matrix(i, j) = sign * det_minor;
        }
    }
    return cofactor_matrix.transpose();
}

} // namespace rational_linalg

#endif // RATIONAL_LINALG_ADJUGATE_FRACTION_HPP
