#pragma once

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>

#include <cstddef>

#include <linalg/integer.hpp>
#include <linalg/matrix_fraction.hpp>

namespace linalg {

/*
 * Owning C++ wrapper around FLINT's row-major arbitrary-precision integer matrix.
 *
 * Entry access returns a non-owning integer reference, so reads and in-place writes use the original FLINT storage without a copy
 * or allocation. The native handle is reserved for the numeric wrapper layer and the specialized fraction-free LDLT kernel.
 */
class matrix_int {
public:
    matrix_int() noexcept { fmpz_mat_init(data_, 0, 0); }

    matrix_int(size_t rows, size_t columns)
    {
        fmpz_mat_init(data_, static_cast<slong>(rows), static_cast<slong>(columns));
    }

    matrix_int(const matrix_int& other) { fmpz_mat_init_set(data_, other.data_); }

    matrix_int(matrix_int&& other) noexcept
    {
        fmpz_mat_init(data_, 0, 0);
        fmpz_mat_swap(data_, other.data_);
    }

    ~matrix_int() noexcept { fmpz_mat_clear(data_); }

    matrix_int& operator=(const matrix_int& other)
    {
        if (this != &other) {
            matrix_int copy(other);
            swap(copy);
        }
        return *this;
    }

    matrix_int& operator=(matrix_int&& other) noexcept
    {
        if (this != &other) fmpz_mat_swap(data_, other.data_);
        return *this;
    }

    size_t rows() const noexcept { return static_cast<size_t>(fmpz_mat_nrows(data_)); }
    size_t cols() const noexcept { return static_cast<size_t>(fmpz_mat_ncols(data_)); }

    integer::reference operator()(size_t row, size_t column) noexcept
    {
        return integer::reference(fmpz_mat_entry(data_, static_cast<slong>(row), static_cast<slong>(column)));
    }

    integer::const_reference operator()(size_t row, size_t column) const noexcept
    {
        return integer::const_reference(fmpz_mat_entry(data_, static_cast<slong>(row), static_cast<slong>(column)));
    }

    void resize(size_t rows, size_t columns)
    {
        if (this->rows() == rows && this->cols() == columns) return;
        matrix_int replacement(rows, columns);
        swap(replacement);
    }

    void negate() noexcept { fmpz_mat_neg(data_, data_); }

    void swap(matrix_int& other) noexcept { fmpz_mat_swap(data_, other.data_); }

    /* Clear one rational matrix's denominators once and retain the minimal positive common denominator separately. */
    void set_from_fraction_matrix(const matrix_frc& source, integer& denominator)
    {
        resize(source.rows(), source.cols());

        fmpq_mat_t rational_matrix;
        fmpq_mat_init(rational_matrix, static_cast<slong>(source.rows()), static_cast<slong>(source.cols()));
        for (size_t row = 0; row < source.rows(); ++row) {
            for (size_t column = 0; column < source.cols(); ++column) {
                fmpq_set(fmpq_mat_entry(rational_matrix, static_cast<slong>(row), static_cast<slong>(column)), source(row, column).data());
            }
        }
        fmpq_mat_get_fmpz_mat_matwise(data_, denominator.native_handle(), rational_matrix);
        fmpq_mat_clear(rational_matrix);
    }

    fmpz_mat_struct* native_handle() noexcept { return data_; }
    const fmpz_mat_struct* native_handle() const noexcept { return data_; }

private:
    fmpz_mat_t data_;
};

inline void swap(matrix_int& left, matrix_int& right) noexcept
{
    left.swap(right);
}

} // namespace linalg
