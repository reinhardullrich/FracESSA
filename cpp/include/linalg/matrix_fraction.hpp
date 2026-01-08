#ifndef RATIONAL_LINALG_MATRIX_FRACTION_HPP
#define RATIONAL_LINALG_MATRIX_FRACTION_HPP

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <linalg/fraction.hpp>
#include <fracessa/bitset64.hpp>

namespace linalg {

class matrix_frc {
public:
    matrix_frc() : rows_(0), cols_(0) {}
    
    // Default constructor: NO initialization (fast - elements must be assigned)
    matrix_frc(size_t rows, size_t cols) 
        : rows_(rows), cols_(cols) {
        data_.resize(rows * cols);
    }

    // In-place initialization methods (no temporaries)
    void set_zero(size_t rows, size_t cols) {
        rows_ = rows;
        cols_ = cols;
        data_.resize(rows * cols);
        for (size_t i = 0; i < rows * cols; ++i) {
            data_[i] = fraction::zero();
        }
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
    
    const std::vector<fraction>& data() const noexcept { return data_; }
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

    /**
     * is_positive_definite - rational LDL^T decomposition.
     * 
     * Accuracy:
     *   - Perfect accuracy with fraction types; no rounding errors.
     *   - Checks if diagonal elements of D are strictly positive.
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

private:
    size_t rows_;
    size_t cols_;
    std::vector<fraction> data_;
};

// Factory functions
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
