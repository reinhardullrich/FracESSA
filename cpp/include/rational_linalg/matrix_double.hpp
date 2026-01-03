#ifndef RATIONAL_LINALG_MATRIX_DOUBLE_HPP
#define RATIONAL_LINALG_MATRIX_DOUBLE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <fracessa/bitset64.hpp>

namespace rational_linalg {

class matrix_double {
public:
    matrix_double() : rows_(0), cols_(0) {}
    matrix_double(size_t rows, size_t cols) 
        : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {}

    static matrix_double identity(size_t n) {
        matrix_double result(n, n);
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = 1.0;
        }
        return result;
    }

    size_t rows() const noexcept { return rows_; }
    size_t cols() const noexcept { return cols_; }
    const std::vector<double>& data() const noexcept { return data_; }
    std::vector<double>& data() noexcept { return data_; }

    double& operator()(size_t i, size_t j) {
        return data_[i * cols_ + j];
    }

    const double& operator()(size_t i, size_t j) const {
        return data_[i * cols_ + j];
    }

    void swap_rows(size_t i, size_t j) {
        if (i == j) return;
        for (size_t k = 0; k < cols_; ++k) {
            std::swap((*this)(i, k), (*this)(j, k));
        }
    }

    matrix_double transpose() const {
        matrix_double result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    matrix_double principal_submatrix(const bitset64& support) const {
        size_t support_size = bs64::count_set_bits(support);
        matrix_double submatrix(support_size, support_size);
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

    double infinity_norm() const {
        double max_norm = 0.0;
        for (size_t i = 0; i < rows_; ++i) {
            double row_sum = 0.0;
            for (size_t j = 0; j < cols_; ++j) {
                row_sum += std::abs((*this)(i, j));
            }
            if (row_sum > max_norm) max_norm = row_sum;
        }
        return max_norm;
    }

    std::string to_log_string() const {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                ss << std::setw(12) << (*this)(i, j) << " ";
            }
            ss << "\n";
        }
        return ss.str();
    }

    // Forward declaration for method defined in separate header
    bool is_positive_definite() const;

private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};

} // namespace rational_linalg

#endif // RATIONAL_LINALG_MATRIX_DOUBLE_HPP
