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

    // Factory function: creates zero-initialized matrix
    static matrix_frc zero(size_t rows, size_t cols) {
        matrix_frc result(rows, cols);
        for (size_t i = 0; i < rows * cols; ++i) {
            result.data_[i] = fraction::zero();
        }
        return result;
    }

    // Factory function: creates identity matrix
    static matrix_frc identity(size_t n) {
        matrix_frc result(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    result(i, j) = fraction::one();
                } else {
                    result(i, j) = fraction::zero();
                }
            }
        }
        return result;
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

    matrix_frc transpose() const {
        matrix_frc result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
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

    bool all_entries_greater_zero() const noexcept {
        for (const auto& val : data_) {
            if (val <= fraction::zero()) return false;
        }
        return true;
    }

    matrix_frc operator*(const matrix_frc& other) const {
        if (cols_ != other.rows_) throw std::runtime_error("Matrix dimensions mismatch");
        matrix_frc result(rows_, other.cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                fraction sum = fraction::zero();
                for (size_t k = 0; k < cols_; ++k) {
                    sum.addmul((*this)(i, k), other(k, j));
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    matrix_frc operator*(const fraction& scalar) const {
        matrix_frc result(rows_, cols_);
        for (size_t i = 0; i < data_.size(); ++i) {
            fraction::mul(result.data_[i], data_[i], scalar);
        }
        return result;
    }

    // Forward declarations for methods defined in separate headers to keep this one lean
    bool is_positive_definite() const;

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
