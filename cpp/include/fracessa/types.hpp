#pragma once

#include <cstddef>
#include <vector>

#include <coposit/integer.hpp>
#include <coposit/matrix_integer.hpp>

namespace fracessa::numeric {

// Zero-cost FracESSA names for the exact types owned by coposit.
using integer = coposit::integer;
using matrix_int = coposit::matrix_integer;

/**
 * Row-major dense double storage for fast candidate search.
 *
 * This is reusable numerical scratch space, not the exact game representation. It deliberately provides only zero-initialized storage,
 * unchecked indexing, and the row count needed by the small reduced systems. Avoiding a generic matrix abstraction keeps these millions of
 * operations cheap.
 */
class matrix_dbl {
public:
    matrix_dbl() : rows_(0), cols_(0) {}
    matrix_dbl(size_t rows, size_t cols)
        : rows_(rows), cols_(cols), data_(rows * cols, 0.0)
    {}

    size_t rows() const noexcept { return rows_; }

    double& operator()(size_t row, size_t column) { return data_[row * cols_ + column]; }
    const double& operator()(size_t row, size_t column) const { return data_[row * cols_ + column]; }

private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};

} // namespace fracessa::numeric
