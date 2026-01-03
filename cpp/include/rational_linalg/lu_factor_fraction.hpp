#ifndef RATIONAL_LINALG_LU_FACTOR_FRACTION_HPP
#define RATIONAL_LINALG_LU_FACTOR_FRACTION_HPP

#include <rational_linalg/matrix_fraction.hpp>
#include <stdexcept>

namespace rational_linalg {

class lu_factor_fraction {
public:
    lu_factor_fraction(const matrix_fraction& A) {
        compute(A);
    }

    void compute(const matrix_fraction& A) {
        const size_t n = A.rows();
        m_n = n;
        m_L = matrix_fraction::identity(n);
        m_U = A;
        m_P = matrix_fraction::identity(n);
        m_swap_count = 0;
        m_is_singular = false;

        if (n == 0) return;

        for (size_t k = 0; k < n - 1; ++k) {
            // Partial Pivoting
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
        if (m_is_singular) return fraction::zero();
        fraction det = fraction::one();
        if (m_swap_count % 2 == 1) det = fraction::neg_one();
        for (size_t i = 0; i < m_n; ++i)
            det *= m_U(i, i);
        return det;
    }

    matrix_fraction inverse() const {
        if (m_is_singular) throw std::runtime_error("Matrix is singular");
        matrix_fraction Inv(m_n, m_n);
        for (size_t col = 0; col < m_n; ++col) {
            matrix_fraction b(m_n, 1);
            b(col, 0) = fraction::one();
            matrix_fraction x = solve(b);
            for (size_t i = 0; i < m_n; ++i) {
                Inv(i, col) = x(i, 0);
            }
        }
        return Inv;
    }

    matrix_fraction solve(const matrix_fraction& b) const {
        
        matrix_fraction bp(m_n, 1);
        for (size_t i = 0; i < m_n; ++i) {
            fraction sum = fraction::zero();
            for (size_t j = 0; j < m_n; ++j) {
                sum.addmul(m_P(i, j), b(j, 0));
            }
            bp(i, 0) = sum;
        }
        
        matrix_fraction y(m_n, 1);
        matrix_fraction x(m_n, 1);

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
    matrix_fraction m_L, m_U, m_P;
    bool m_is_singular = false;
    int m_swap_count = 0;
};

} // namespace rational_linalg

#endif // RATIONAL_LINALG_LU_FACTOR_FRACTION_HPP
