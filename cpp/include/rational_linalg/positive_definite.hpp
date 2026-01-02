#pragma once

#include <rational_linalg/matrix.hpp>
#include <rational_linalg/fraction.hpp>
#include <rational_linalg/constants.hpp>
#include <type_traits>
#include <cmath>
#include <limits>

namespace rational_linalg {

// Member function implementation for is_positive_definite
// Uses LDLT decomposition for fraction type and Cholesky decomposition for double type
template<typename T>
inline bool Matrix<T>::is_positive_definite() const {
    if constexpr (std::is_same_v<T, fraction>) {
        // LDLT decomposition for fractions
        size_t n = rows_;
        Matrix<T> D(n, n);
        Matrix<T> L(n, n);
        
        // Initialize L as identity
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                L(i, j) = (i == j) ? rational_linalg::one<T>() : rational_linalg::zero<T>();
            }
        }

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                T aSum = rational_linalg::zero<T>();
                for (size_t k = 0; k < j; ++k) {
                    aSum += L(i, k) * L(j, k) * D(k, k);
                }
                L(i, j) = (rational_linalg::one<T>() / D(j, j)) * ((*this)(i, j) - aSum);
            }
            T bSum = rational_linalg::zero<T>();
            for (size_t k = 0; k < i; ++k) {
                bSum += L(i, k) * L(i, k) * D(k, k);
            }
            D(i, i) = (*this)(i, i) - bSum;
            if (D(i, i) <= rational_linalg::zero<T>()) {
                return false;
            }
        }
        return true;
    } else if constexpr (std::is_same_v<T, double>) {
        // Cholesky decomposition for double matrices
        const size_t n = rows_;
        
        // Matrix must be square
        if (n != cols_) {
            return false;
        }
        
        // Initialize L as zero matrix
        Matrix<double> L = Matrix<double>::zero(n, n);
        
        // Tolerance for positive definiteness check
        // Use machine epsilon scaled by infinity norm for numerical stability
        const double tolerance = std::numeric_limits<double>::epsilon() * infinity_norm();
        
        for (size_t i = 0; i < n; ++i) {
            // Compute off-diagonal elements: L[i,j] for j < i
            for (size_t j = 0; j < i; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < j; ++k) {
                    sum += L(i, k) * L(j, k);
                }
                L(i, j) = (1.0 / L(j, j)) * ((*this)(i, j) - sum);
            }
            
            // Compute diagonal element: L[i,i]
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L(i, k) * L(i, k);
            }
            double diagonal_value = (*this)(i, i) - sum;
            
            // Check if matrix is positive definite
            // Diagonal value must be positive (within tolerance)
            if (diagonal_value <= tolerance) {
                return false;
            }
            
            L(i, i) = std::sqrt(diagonal_value);
        }
        
        return true;
    }
} 
}// namespace rational_linalg

