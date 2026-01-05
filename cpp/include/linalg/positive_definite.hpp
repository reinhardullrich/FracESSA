#ifndef RATIONAL_LINALG_POSITIVE_DEFINITE_HPP
#define RATIONAL_LINALG_POSITIVE_DEFINITE_HPP

#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_double.hpp>
#include <limits>
#include <cmath>
#include <vector>

namespace linalg {

// ============================================================================
// Double Implementation
// ============================================================================

/**
 * is_positive_definite - Cholesky decomposition for double matrices.
 * 
 * Stability / Accuracy:
 *   - Uses machine epsilon scaled by infinity norm as a tolerance for 
 *     better numerical stability in edge cases.
 */
inline bool matrix_dbl::is_positive_definite() const {
    const size_t n = rows_;
    if (n == 0) return true;
    if (n != cols_) return false;
    
    matrix_dbl L(n, n);
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

// ============================================================================
// Fraction Implementation
// ============================================================================

/**
 * is_positive_definite - rational LDL^T decomposition.
 * 
 * Accuracy:
 *   - Perfect accuracy with fraction types; no rounding errors.
 *   - Checks if diagonal elements of D are strictly positive.
 */
inline bool matrix_frc::is_positive_definite() const {
    const size_t n = rows_;
    if (n == 0) return true;
    
    std::vector<fraction> d(n);
    matrix_frc L = matrix_frc::zero(n, n);
    
    for (size_t i = 0; i < n; ++i) {
        L(i, i) = fraction::one();
    }

    fraction aSum, bSum, term;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                aSum = fraction::zero();
                for (size_t k = 0; k < j; ++k) {
                    // term = L(j, k) * d[k]
                    fraction::mul(term, L(j, k), d[k]);
                    // aSum += L(i, k) * term
                    aSum.addmul(L(i, k), term);
                }
                // L(i, j) = A(i, j) - aSum
                fraction::sub(L(i, j), (*this)(i, j), aSum);
                // L(i, j) /= d[j]
                L(i, j).div_inplace(d[j]);
            }
            
            bSum = fraction::zero();
        for (size_t k = 0; k < i; ++k) {
            // term = L(i, k) * d[k]
            fraction::mul(term, L(i, k), d[k]);
            // bSum += L(i, k) * term
            bSum.addmul(L(i, k), term);
        }
        
        // d[i] = A(i, i) - bSum
        fraction::sub(d[i], (*this)(i, i), bSum);
        if (d[i].sgn() <= 0) {
            return false;
        }
    }
    return true;
}

} // namespace linalg

#endif // RATIONAL_LINALG_POSITIVE_DEFINITE_HPP
