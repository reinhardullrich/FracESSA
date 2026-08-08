#ifndef RATIONAL_LINALG_MATRIX_FRACTION_HPP
#define RATIONAL_LINALG_MATRIX_FRACTION_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <linalg/fraction.hpp>

namespace linalg {

/*
 * Dense exact-rational storage for parsed games and candidate output.
 *
 * Design intent:
 * Exact candidate and stability calculations use integer matrices; this type
 * remains only at the rational input/output boundary.
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

    size_t rows() const noexcept { return rows_; }
    size_t cols() const noexcept { return cols_; }

    fraction& operator()(size_t i, size_t j) {
        return data_[i * cols_ + j];
    }

    const fraction& operator()(size_t i, size_t j) const {
        return data_[i * cols_ + j];
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
