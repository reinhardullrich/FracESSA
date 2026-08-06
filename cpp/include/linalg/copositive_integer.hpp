#ifndef LINALG_COPOSITIVITY_INTEGER_HPP
#define LINALG_COPOSITIVITY_INTEGER_HPP

#include <linalg/fraction_free_ldlt.hpp>
#include <linalg/matrix_integer.hpp>
#include <fracessa/bitset64.hpp>
#include <array>
#include <cstddef>

namespace linalg {

/*
 * Facts obtained from one exact sign scan of a symmetric matrix. The negative
 * neighbor masks are the graph needed by the later component reduction. The
 * row masks and negative-part row sums answer the remaining sign questions
 * without another matrix pass.
 */
struct CopositivitySignScan {
    // Caller contract: test all_diagonal_positive before reading any other field. A false value means the scan returned before the
    // off-diagonal pass, so every other field is intentionally incomplete.
    std::array<bitset64, bs64::kMaxBitsetDimension> negative_neighbors{};
    bitset64 rows_with_negative_off_diagonal = 0;
    bitset64 rows_with_positive_off_diagonal = 0;
    integer all_ones_quadratic_value;
    bool all_diagonal_positive = true;
    bool all_negative_part_row_sums_positive = true;
};

/*
 * Scan the integer-scaled reduced B matrix once. The quadratic sum and every
 * sign decision remain exact even when entries exceed a machine word.
 */
inline CopositivitySignScan scan_copositivity_signs(const matrix_int& A) {
    CopositivitySignScan result;
    std::array<integer, bs64::kMaxBitsetDimension> negative_part_row_sums;
    const size_t n = A.rows();

    // This must remain the first pass: a failed diagonal returns an incomplete result whose other fields must not be consumed.
    for (size_t i = 0; i < n; ++i) {
        const integer::const_reference diagonal = A(i, i);
        if (diagonal.sign() <= 0) {
            result.all_diagonal_positive = false;
            return result;
        }
        negative_part_row_sums[i] = diagonal;
        result.all_ones_quadratic_value += diagonal;
    }

    // Symmetry lets one triangular scan build both graph directions and count
    // every off-diagonal contribution to 1^T A 1 exactly twice.
    for (size_t i = 0; i < n; ++i) {
        const bitset64 bit_i = bs64::single_bit_at_pos(i);
        for (size_t j = i + 1; j < n; ++j) {
            const integer::const_reference entry = A(i, j);
            const bitset64 bit_j = bs64::single_bit_at_pos(j);
            const bitset64 endpoints = bit_i | bit_j;
            const int sign = entry.sign();

            if (sign < 0) {
                result.negative_neighbors[i] |= bit_j;
                result.negative_neighbors[j] |= bit_i;
                result.rows_with_negative_off_diagonal |= endpoints;
                negative_part_row_sums[i] += entry;
                negative_part_row_sums[j] += entry;
            } else if (sign > 0) {
                result.rows_with_positive_off_diagonal |= endpoints;
            }

            result.all_ones_quadratic_value += entry;
            result.all_ones_quadratic_value += entry;
        }

        if (negative_part_row_sums[i].sign() <= 0) result.all_negative_part_row_sums_positive = false;
    }

    return result;
}

// Exact low-dimensional criteria from Bomze (1992), Theorem 3.4. Only one triangle is needed because B is symmetric.
inline bool is_strictly_copositive_1x1(integer::const_reference b11) noexcept {
    return b11.sign() > 0;
}

inline bool is_strictly_copositive_2x2(integer::const_reference b11, integer::const_reference b12,
                                      integer::const_reference b22) noexcept {
    if (b11.sign() <= 0 || b22.sign() <= 0) return false;
    if (b12.sign() >= 0) return true;

    integer determinant;
    determinant.set_product(b11, b22);
    determinant.submul(b12, b12);
    return determinant.sign() > 0;
}

inline bool is_strictly_copositive_3x3(integer::const_reference b11, integer::const_reference b12,
                                       integer::const_reference b13, integer::const_reference b22,
                                       integer::const_reference b23, integer::const_reference b33) noexcept {
    if (b11.sign() <= 0 || b22.sign() <= 0 || b33.sign() <= 0) return false;

    // Every 2x2 principal submatrix must first be strictly copositive.
    integer work;
    if (b12.sign() < 0) {
        work.set_product(b11, b22);
        work.submul(b12, b12);
        if (work.sign() <= 0) return false;
    }
    if (b13.sign() < 0) {
        work.set_product(b11, b33);
        work.submul(b13, b13);
        if (work.sign() <= 0) return false;
    }
    if (b23.sign() < 0) {
        work.set_product(b22, b33);
        work.submul(b23, b23);
        if (work.sign() <= 0) return false;
    }

    // det(B) = b11*b22*b33 + 2*b12*b13*b23 - b11*b23^2 - b22*b13^2 - b33*b12^2.
    integer determinant;
    determinant.set_product(b11, b22);
    fmpz_mul(determinant.native_handle(), determinant.native_handle(), b33.native_handle());
    work.set_product(b12, b13);
    fmpz_mul(work.native_handle(), work.native_handle(), b23.native_handle());
    work.multiply(2);
    determinant += work;
    work.set_product(b23, b23);
    determinant.submul(b11, work);
    work.set_product(b13, b13);
    determinant.submul(b22, work);
    work.set_product(b12, b12);
    determinant.submul(b33, work);

    if (determinant.sign() > 0) return true;

    // Hadeler rejects exactly when all six distinct entries of the symmetric adjugate are positive.
    work.set_product(b22, b33);
    work.submul(b23, b23);
    if (work.sign() <= 0) return true;
    work.set_product(b11, b33);
    work.submul(b13, b13);
    if (work.sign() <= 0) return true;
    work.set_product(b11, b22);
    work.submul(b12, b12);
    if (work.sign() <= 0) return true;
    work.set_product(b13, b23);
    work.submul(b12, b33);
    if (work.sign() <= 0) return true;
    work.set_product(b12, b23);
    work.submul(b13, b22);
    if (work.sign() <= 0) return true;
    work.set_product(b12, b13);
    work.submul(b11, b23);
    return work.sign() <= 0;
}

/*
 * Exact strict-copositivity test for a symmetric integer matrix A:
 *
 *     z^T A z > 0 for every nonzero z >= 0.
 *
 * FracESSA reaches this code only after the simpler pure-ESS and exact
 * positive-definiteness cases have failed. It checks principal submatrices in
 * increasing order. Therefore, when a matrix of size k is reached, all of its
 * proper principal submatrices have already passed. Hadeler (1983), Theorem 3,
 * then decides the remaining case from the determinant and adjugate signs.
 */
class CopositivityChecker {
private:
    // Compute adj(A[subset,subset]) from cofactors when the principal matrix is singular and its inverse is
    // unavailable.
    matrix_int adjugate(const matrix_int& A, const uint8_t* subset_indices, size_t n) {
        matrix_int adj(n, n);
        matrix_int minor(n - 1, n - 1);

        // A symmetric matrix has a symmetric adjugate, so calculate one triangle.
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                // Extract the minor obtained by deleting row i and column j.
                size_t minor_row = 0;
                for (size_t row = 0; row < n; ++row) {
                    if (row == i) continue;
                    size_t minor_col = 0;
                    for (size_t col = 0; col < n; ++col) {
                        if (col == j) continue;
                        minor(minor_row, minor_col) = A(subset_indices[row], subset_indices[col]);
                        ++minor_col;
                    }
                    ++minor_row;
                }

                // Off-diagonal cofactor minors are not symmetric, so the symmetric LDL^T factorization cannot
                // calculate them. FLINT's exact integer determinant replaces the former rational LU call without
                // changing the cofactor formula.
                integer::reference cofactor = adj(j, i);
                fmpz_mat_det(cofactor.native_handle(), minor.native_handle());
                if (((i + j) & 1) != 0) cofactor.negate();

                // adj(A) is the transposed cofactor matrix; symmetry makes both
                // mirrored entries equal.
                if (i != j) adj(i, j) = cofactor;
            }
        }
        return adj;
    }

    bool all_entries_greater_zero(const matrix_int& A) const noexcept {
        const size_t n = A.rows();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                // Only one triangle is needed because every input here is symmetric.
                if (A(i, j).sign() <= 0) return false;
            }
        }
        return true;
    }

    bool is_copositive_hadeler(const matrix_int& A, bitset64 mask, size_t current_dim) {
        switch (current_dim) {
            case 1: {
                const size_t i = bs64::find_pos_first_set_bit(mask);
                return is_strictly_copositive_1x1(A(i, i));
            }
            case 2: {
                uint8_t indices[bs64::kMaxBitsetDimension];
                bs64::extract_set_indices(mask, A.rows(), indices);
                const size_t i = indices[0];
                const size_t j = indices[1];
                return is_strictly_copositive_2x2(A(i, i), A(i, j), A(j, j));
            }
            case 3: {
                uint8_t indices[bs64::kMaxBitsetDimension];
                bs64::extract_set_indices(mask, A.rows(), indices);
                const size_t i = indices[0];
                const size_t j = indices[1];
                const size_t k = indices[2];
                return is_strictly_copositive_3x3(A(i, i), A(i, j), A(i, k), A(j, j), A(j, k), A(k, k));
            }
            default:
                break;
        }

        uint8_t subset_indices[bs64::kMaxBitsetDimension];
        bs64::extract_set_indices(mask, A.rows(), subset_indices);
        matrix_int subMat(current_dim, current_dim);
        for (size_t row = 0; row < current_dim; ++row) {
            for (size_t column = 0; column < current_dim; ++column) {
                subMat(row, column) = A(subset_indices[row], subset_indices[column]);
            }
        }

        const bool nonsingular = factorization_.factorize_inplace(subMat) != 0;
        const int determinant_sign = factorization_.determinant().sign();

        if (determinant_sign <= 0) {
            /*
             * With all proper principal submatrices already passed, Hadeler's
             * rejecting pattern is det(A) <= 0 and adj(A) > 0 entrywise. In the
             * singular case, a positive adjugate gives a positive null direction,
             * so the quadratic form is not strictly positive there.
             */
            matrix_int adj;
            if (!nonsingular) {
                // No inverse exists, so use the cofactor definition.
                adj = adjugate(A, subset_indices, current_dim);
            } else {
                // solve_inplace returns |det(A)|*A^-1. This branch has det(A)<0, so its negation is adj(A).
                adj.set_identity(current_dim);
                integer denominator;
                factorization_.solve_inplace(adj, denominator, subMat);
                adj.negate();
            }

            if (all_entries_greater_zero(adj)) {
                return false;
            }
        }

        return true;
    }

public:
    explicit CopositivityChecker(size_t maximum_dimension)
        : factorization_(maximum_dimension)
    {
    }

    bool is_positive_definite(const matrix_int& A) {
        matrix_int factored(A);
        return factorization_.factorize_inplace(factored) != 0 && factorization_.is_positive_definite();
    }

    bool is_strictly_copositive(const matrix_int& A) {
        size_t n = A.rows();
        const bitset64 limit = bs64::two_to_the_power_of(n);

        // Cardinality order establishes Hadeler's proper-submatrix precondition.
        // Gosper's algorithm streams each layer and permits immediate rejection.
        for (size_t subset_size = 1; subset_size <= n; ++subset_size) {
            bitset64 subset = bs64::set_all_n_bits(subset_size);
            while (subset < limit) {
                if (!is_copositive_hadeler(A, subset, subset_size)) {
                    return false;
                }

                subset = bs64::next_same_cardinality(subset);
            }
        }

        // In particular, the final full-size subset passed.
        return true;
    }

private:
    fraction_free_ldlt_factorization factorization_;
};

inline bool is_strictly_copositive(const matrix_int& A) {
    CopositivityChecker checker(A.rows());
    return checker.is_strictly_copositive(A);
}

} // namespace linalg

#endif // LINALG_COPOSITIVITY_INTEGER_HPP
