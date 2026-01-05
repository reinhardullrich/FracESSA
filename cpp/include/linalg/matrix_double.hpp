#ifndef RATIONAL_LINALG_MATRIX_DOUBLE_HPP
#define RATIONAL_LINALG_MATRIX_DOUBLE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <fracessa/bitset64.hpp>

namespace linalg {

class matrix_dbl {
public:
    matrix_dbl() : rows_(0), cols_(0) {}
    matrix_dbl(size_t rows, size_t cols) 
        : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {}

    static matrix_dbl identity(size_t n) {
        matrix_dbl result(n, n);
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

    matrix_dbl transpose() const {
        matrix_dbl result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
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

} // namespace linalg

#endif // RATIONAL_LINALG_MATRIX_DOUBLE_HPP
