#ifndef RATIONAL_LINALG_MATRIX_FRACTION_HPP
#define RATIONAL_LINALG_MATRIX_FRACTION_HPP

#include <cstddef>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <linalg/fraction.hpp>

namespace linalg {

/*
 * Dense exact-rational matrix for FracESSA kernels.
 *
 * Design intent:
 * - Keep storage contiguous for cache-friendly scans in elimination/factorization.
 * - Avoid generic linear-algebra abstraction overhead; only operations needed by
 *   the ESS pipeline are implemented.
 */
class matrix_frc {
public:
    matrix_frc() : rows_(0), cols_(0) {}
    
    // Constructor resizes storage; elements become default-constructed fractions.
    // Callers in hot paths usually overwrite all entries immediately.
    matrix_frc(size_t rows, size_t cols) 
        : rows_(rows), cols_(cols) {
        data_.resize(rows * cols);
    }

    void set_identity(size_t n) {
        rows_ = n;
        cols_ = n;
        data_.resize(n * n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    data_[i * n + j] = fraction::one();
                } else {
                    data_[i * n + j] = fraction::zero();
                }
            }
        }
    }

    size_t rows() const noexcept { return rows_; }
    size_t cols() const noexcept { return cols_; }
    
    std::vector<fraction>& data() noexcept { return data_; }

    fraction& operator()(size_t i, size_t j) {
        return data_[i * cols_ + j];
    }

    const fraction& operator()(size_t i, size_t j) const {
        return data_[i * cols_ + j];
    }

    void swap_rows(size_t i, size_t j) {
        for (size_t k = 0; k < cols_; ++k) {
            std::swap((*this)(i, k), (*this)(j, k));
        }
    }


    std::string to_log_string() const {
        std::stringstream ss;
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                ss << std::setw(12) << (*this)(i, j).to_string() << " ";
            }
            ss << "\n";
        }
        return ss.str();
    }


    matrix_frc operator*(const fraction& scalar) const {
        matrix_frc result(rows_, cols_);
        for (size_t i = 0; i < data_.size(); ++i) {
            fraction::mul(result.data_[i], data_[i], scalar);
        }
        return result;
    }

    /*
     * Test a symmetric matrix by the exact factorization A = L*D*L^T, where L
     * has unit diagonal and D is diagonal. A is positive definite exactly when
     * every entry of D is positive. Because all entries are rational, these are
     * exact sign decisions: no epsilon or floating-point tolerance is involved.
     * This is the fast sufficient stability test for the Bee matrix.
     */
    bool is_positive_definite() const {
        const size_t n = rows_;
        if (n == 0) return true;
        
        std::vector<fraction> d(n);
        
        matrix_frc L;
        L.set_identity(n);

        fraction aSum, bSum, term;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                aSum = fraction::zero();
                for (size_t k = 0; k < j; ++k) {
                    // Subtract the part already explained by earlier columns of L.
                    fraction::mul(term, L(j, k), d[k]);
                    aSum.addmul(L(i, k), term);
                }
                fraction::sub(L(i, j), (*this)(i, j), aSum);
                L(i, j).div_inplace(d[j]);
            }
            
            bSum = fraction::zero();
            for (size_t k = 0; k < i; ++k) {
                fraction::mul(term, L(i, k), d[k]);
                bSum.addmul(L(i, k), term);
            }
            
            // This is the next diagonal entry of D; its sign decides this step.
            fraction::sub(d[i], (*this)(i, i), bSum);
            if (d[i].sgn() <= 0) {
                return false;
            }
        }
        return true;
    }

private:
    size_t rows_;
    size_t cols_;
    std::vector<fraction> data_;
};

/*
 * Expand the compact circular-symmetric format. `half_row` gives the payoff by
 * circular distance from strategy 0; reflection supplies the opposite side,
 * and every later row is a cyclic shift of the first. The diagonal is zero by
 * definition of this input format.
 */
inline matrix_frc create_circular_symmetric(size_t n, const std::vector<fraction>& half_row) {
    matrix_frc result(n, n);
    std::vector<fraction> first_row(n);
    first_row[0] = fraction::zero();
    
    if (n % 2 == 0) {
        first_row[n / 2] = half_row[half_row.size() - 1];
        for (size_t i = 0; i < n / 2 - 1; ++i) {
            first_row[i + 1] = half_row[i];
            first_row[n - i - 1] = half_row[i];
        }
    } else {
        for (size_t i = 0; i < n / 2; ++i) {
            first_row[i + 1] = half_row[i];
            first_row[n - i - 1] = half_row[i];
        }
    }
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result(i, j) = first_row[(j - i + n) % n];
        }
    }
    return result;
}

// Expand row-major upper-triangle input and mirror it across the diagonal.
inline matrix_frc create_symmetric(size_t n, const std::vector<fraction>& upper_triangular) {
    matrix_frc result(n, n);
    size_t idx = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            fraction val = upper_triangular[idx];
            result(i, j) = val;
            result(j, i) = val;
            ++idx;
        }
    }
    return result;
}

} // namespace linalg

#endif // RATIONAL_LINALG_MATRIX_FRACTION_HPP
