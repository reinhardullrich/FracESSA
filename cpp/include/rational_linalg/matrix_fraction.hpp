#ifndef RATIONAL_LINALG_MATRIX_FRACTION_HPP
#define RATIONAL_LINALG_MATRIX_FRACTION_HPP

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <rational_linalg/fraction.hpp>
#include <fracessa/bitset64.hpp>

namespace rational_linalg {

class matrix_fraction {
public:
    matrix_fraction() : rows_(0), cols_(0) {}
    matrix_fraction(size_t rows, size_t cols) 
        : rows_(rows), cols_(cols), data_(rows * cols, fraction::zero()) {}

    static matrix_fraction identity(size_t n) {
        matrix_fraction result(n, n);
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = fraction::one();
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

    matrix_fraction transpose() const {
        matrix_fraction result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    matrix_fraction principal_submatrix(const bitset64& support) const {
        size_t support_size = bs64::count_set_bits(support);
        matrix_fraction submatrix(support_size, support_size);
        size_t row = 0;
        for (size_t i = bs64::find_pos_first_set_bit(support); i < 64; i = bs64::find_pos_next_set_bit(support, i)) {
            size_t col = 0;
            for (size_t j = bs64::find_pos_first_set_bit(support); j < 64; j = bs64::find_pos_next_set_bit(support, j)) {
                submatrix(row, col) = (*this)(i, j);
                ++col;
            }
            ++row;
        }
        return submatrix;
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

    matrix_fraction operator*(const matrix_fraction& other) const {
        if (cols_ != other.rows_) throw std::runtime_error("Matrix dimensions mismatch");
        matrix_fraction result(rows_, other.cols_);
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

    matrix_fraction operator*(const fraction& scalar) const {
        matrix_fraction result(rows_, cols_);
        for (size_t i = 0; i < data_.size(); ++i) {
            fraction::mul(result.data_[i], data_[i], scalar);
        }
        return result;
    }

    // Forward declarations for methods defined in separate headers to keep this one lean
    bool is_positive_definite() const;
    std::vector<double> to_double_vec() const;

private:
    size_t rows_;
    size_t cols_;
    std::vector<fraction> data_;
};

// Factory functions
inline matrix_fraction create_circular_symmetric(size_t n, const std::vector<fraction>& half_row) {
    matrix_fraction result(n, n);
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

inline matrix_fraction create_symmetric(size_t n, const std::vector<fraction>& upper_triangular) {
    matrix_fraction result(n, n);
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

} // namespace rational_linalg

#endif // RATIONAL_LINALG_MATRIX_FRACTION_HPP
