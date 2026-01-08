#ifndef RATIONAL_LINALG_MATRIX_DOUBLE_HPP
#define RATIONAL_LINALG_MATRIX_DOUBLE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <fracessa/bitset64.hpp>

namespace linalg {

class matrix_dbl {
public:
    matrix_dbl() : rows_(0), cols_(0) {}
    matrix_dbl(size_t rows, size_t cols) 
        : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {}


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

    /**
     * is_positive_definite - Cholesky decomposition for double matrices.
     * 
     * Stability / Accuracy:
     *   - Uses machine epsilon scaled by infinity norm as a tolerance for 
     *     better numerical stability in edge cases.
     */
    bool is_positive_definite() const {
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

private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};

} // namespace linalg

#endif // RATIONAL_LINALG_MATRIX_DOUBLE_HPP
