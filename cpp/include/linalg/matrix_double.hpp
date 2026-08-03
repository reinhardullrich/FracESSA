#ifndef RATIONAL_LINALG_MATRIX_DOUBLE_HPP
#define RATIONAL_LINALG_MATRIX_DOUBLE_HPP

#include <cstddef>
#include <vector>

namespace linalg {

/*
 * Row-major dense double storage for fast candidate search.
 *
 * This is reusable numerical scratch space, not the exact game representation.
 * It deliberately provides only zero-initialized storage, unchecked indexing,
 * and the row count needed by the small reduced systems. Avoiding a generic
 * matrix abstraction keeps these millions of operations cheap.
 */
class matrix_dbl {
public:
    matrix_dbl() : rows_(0), cols_(0) {}
    matrix_dbl(size_t rows, size_t cols) 
        : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {}


    size_t rows() const noexcept { return rows_; }

    double& operator()(size_t i, size_t j) {
        return data_[i * cols_ + j];
    }

    const double& operator()(size_t i, size_t j) const {
        return data_[i * cols_ + j];
    }

private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};

} // namespace linalg

#endif // RATIONAL_LINALG_MATRIX_DOUBLE_HPP
