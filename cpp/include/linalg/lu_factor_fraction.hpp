#ifndef RATIONAL_LINALG_LU_FACTOR_FRACTION_HPP
#define RATIONAL_LINALG_LU_FACTOR_FRACTION_HPP

#include <linalg/matrix_fraction.hpp>
#include <stdexcept>

namespace linalg {

/*
 * Exact LU factorization (with partial pivoting) for rational matrices.
 *
 * This module is used by copositivity checks for determinant, inverse and
 * solve operations on principal submatrices, where exact signs are critical.
 */
/**
 * LU_Factorization (Standard LU factorization with partial pivoting)
 *
 * Performance:
 *   - Complexity: O(n^3), similar to Bareiss LU.
 *   - Efficiency: In rational arithmetic, standard LU is often faster than Bareiss 
 *     because intermediate coefficients grow less rapidly when GCD-canonicalized.
 *
 * Stability / Accuracy:
 *   - Perfect for rational types; no rounding errors.
 *   - Detects singular matrices exactly (pivot = 0).
 */
class LU_Factorization {
public:
    /**
     * Constructor performs a DEEP COPY of A into internal storage m_U.
     * Original matrix A is NOT modified.
     */
    explicit LU_Factorization(const matrix_frc& A) {
        compute(A);
    }

    void compute(const matrix_frc& A) {
        const size_t n = A.rows();
        m_n = n;
        m_L.set_identity(n);
        m_U = A;
        m_P.set_identity(n);
        m_swap_count = 0;
        m_is_singular = false;

        if (n == 0) return;

        // Standard exact P*A = L*U with partial pivoting.
        for (size_t k = 0; k < n - 1; ++k) {
            // Partial pivoting by largest absolute entry in current column.
            size_t max_row = k;
            fraction max_val = m_U(k, k);
            if (max_val < fraction::zero()) max_val = -max_val;
            
            for (size_t i = k + 1; i < n; ++i) {
                fraction val = m_U(i, k);
                if (val < fraction::zero()) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            if (max_row != k) {
                m_U.swap_rows(k, max_row);
                m_P.swap_rows(k, max_row);
                if (k > 0) {
                    for (size_t j = 0; j < k; ++j) {
                        std::swap(m_L(k, j), m_L(max_row, j));
                    }
                }
                m_swap_count++;
            }
            
            const fraction pivot = m_U(k, k);
            // Exact singularity detection (no tolerance required in rationals).
            if (pivot == fraction::zero()) {
                m_is_singular = true;
                return;
            }

            for (size_t i = k + 1; i < n; ++i) {
                m_L(i, k) = m_U(i, k) / pivot;
                for (size_t j = k + 1; j < n; ++j) {
                    m_U(i, j).submul(m_L(i, k), m_U(k, j));
                }
                m_U(i, k) = fraction::zero();
            }
        }

        m_is_singular = (m_U(n - 1, n - 1) == fraction::zero());
    }

    fraction determinant() const {
        // det(P*A) = det(L)*det(U), with det(L)=1 for unit-lower-triangular L.
        if (m_is_singular) return fraction::zero();
        fraction det = fraction::one();
        if (m_swap_count % 2 == 1) det = fraction::neg_one();
        for (size_t i = 0; i < m_n; ++i)
            det *= m_U(i, i);
        return det;
    }

    matrix_frc inverse() const {
        if (m_is_singular) throw std::runtime_error("Matrix is singular");
        matrix_frc Inv(m_n, m_n);
        for (size_t col = 0; col < m_n; ++col) {
            matrix_frc b(m_n, 1);
            b(col, 0) = fraction::one();
            matrix_frc x = solve(b);
            for (size_t i = 0; i < m_n; ++i) {
                Inv(i, col) = x(i, 0);
            }
        }
        return Inv;
    }

    matrix_frc solve(const matrix_frc& b) const {
        // Apply permutation P first, then forward/back substitution.
        matrix_frc bp(m_n, 1);
        for (size_t i = 0; i < m_n; ++i) {
            fraction sum = fraction::zero();
            for (size_t j = 0; j < m_n; ++j) {
                sum.addmul(m_P(i, j), b(j, 0));
            }
            bp(i, 0) = sum;
        }
        
        matrix_frc y(m_n, 1);
        matrix_frc x(m_n, 1);

        for (size_t i = 0; i < m_n; ++i) {
            fraction sum = bp(i, 0);
            for (size_t j = 0; j < i; ++j)
                sum.submul(m_L(i, j), y(j, 0));
            y(i, 0) = sum;
        }

        for (size_t i = m_n; i-- > 0; ) {
            fraction sum = y(i, 0);
            for (size_t j = i + 1; j < m_n; ++j)
                sum.submul(m_U(i, j), x(j, 0));
            x(i, 0) = sum / m_U(i, i);
        }

        return x;
    }

    bool isSingular() const { return m_is_singular; }

private:
    size_t m_n;
    matrix_frc m_L, m_U, m_P;
    bool m_is_singular = false;
    int m_swap_count = 0;
};

} // namespace linalg

#endif // RATIONAL_LINALG_LU_FACTOR_FRACTION_HPP
