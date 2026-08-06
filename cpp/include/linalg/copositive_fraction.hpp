#ifndef RATIONAL_LINALG_COPOSITIVITY_FRACTION_HPP
#define RATIONAL_LINALG_COPOSITIVITY_FRACTION_HPP

#include <linalg/matrix_fraction.hpp>
#include <linalg/lu_factor_fraction.hpp>
#include <fracessa/bitset64.hpp>
#include <array>
#include <cassert>
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
 * Scan the denominator-one scaled reduced B matrix. Its entries are stored as
 * fractions only for the existing matrix API, so the quadratic sum adds their
 * integer numerators directly.
 */
inline CopositivitySignScan scan_copositivity_signs(const matrix_frc& A) {
    CopositivitySignScan result;
    std::array<integer, bs64::kMaxBitsetDimension> negative_part_row_sums;
    const size_t n = A.rows();

    // This must remain the first pass: a failed diagonal returns an incomplete result whose other fields must not be consumed.
    for (size_t i = 0; i < n; ++i) {
        const fraction& diagonal = A(i, i);
        assert(fmpz_is_one(fmpq_denref(diagonal.data())));
        if (diagonal.sgn() <= 0) {
            result.all_diagonal_positive = false;
            return result;
        }
        const integer::const_reference numerator(fmpq_numref(diagonal.data()));
        negative_part_row_sums[i] = numerator;
        result.all_ones_quadratic_value += numerator;
    }

    // Symmetry lets one triangular scan build both graph directions and count
    // every off-diagonal contribution to 1^T A 1 exactly twice.
    for (size_t i = 0; i < n; ++i) {
        const bitset64 bit_i = bs64::single_bit_at_pos(i);
        for (size_t j = i + 1; j < n; ++j) {
            const fraction& entry = A(i, j);
            assert(fmpz_is_one(fmpq_denref(entry.data())));
            const bitset64 bit_j = bs64::single_bit_at_pos(j);
            const bitset64 endpoints = bit_i | bit_j;
            const int sign = entry.sgn();
            const integer::const_reference numerator(fmpq_numref(entry.data()));

            if (sign < 0) {
                result.negative_neighbors[i] |= bit_j;
                result.negative_neighbors[j] |= bit_i;
                result.rows_with_negative_off_diagonal |= endpoints;
                negative_part_row_sums[i] += numerator;
                negative_part_row_sums[j] += numerator;
            } else if (sign > 0) {
                result.rows_with_positive_off_diagonal |= endpoints;
            }

            result.all_ones_quadratic_value += numerator;
            result.all_ones_quadratic_value += numerator;
        }

        if (negative_part_row_sums[i].sign() <= 0) result.all_negative_part_row_sums_positive = false;
    }

    return result;
}

// Exact low-dimensional criteria from Bomze (1992), Theorem 3.4. Only one triangle is needed because B is symmetric.
inline bool is_strictly_copositive_1x1(const fraction& b11) noexcept {
    return b11.sgn() > 0;
}

inline bool is_strictly_copositive_2x2(const fraction& b11, const fraction& b12, const fraction& b22) noexcept {
    if (b11.sgn() <= 0 || b22.sgn() <= 0) return false;
    if (b12.sgn() >= 0) return true;

    fraction determinant;
    fraction::mul(determinant, b11, b22);
    determinant.submul(b12, b12);
    return determinant.sgn() > 0;
}

inline bool is_strictly_copositive_3x3(const fraction& b11, const fraction& b12, const fraction& b13,
                                       const fraction& b22, const fraction& b23, const fraction& b33) noexcept {
    if (b11.sgn() <= 0 || b22.sgn() <= 0 || b33.sgn() <= 0) return false;

    // Every 2x2 principal submatrix must first be strictly copositive.
    fraction work;
    if (b12.sgn() < 0) {
        fraction::mul(work, b11, b22);
        work.submul(b12, b12);
        if (work.sgn() <= 0) return false;
    }
    if (b13.sgn() < 0) {
        fraction::mul(work, b11, b33);
        work.submul(b13, b13);
        if (work.sgn() <= 0) return false;
    }
    if (b23.sgn() < 0) {
        fraction::mul(work, b22, b33);
        work.submul(b23, b23);
        if (work.sgn() <= 0) return false;
    }

    // det(B) = b11*b22*b33 + 2*b12*b13*b23 - b11*b23^2 - b22*b13^2 - b33*b12^2.
    fraction determinant;
    fraction::mul(determinant, b11, b22);
    determinant *= b33;
    fraction::mul(work, b12, b13);
    work *= b23;
    determinant.addmul(work, fraction::two());
    fraction::mul(work, b23, b23);
    determinant.submul(b11, work);
    fraction::mul(work, b13, b13);
    determinant.submul(b22, work);
    fraction::mul(work, b12, b12);
    determinant.submul(b33, work);

    if (determinant.sgn() > 0) return true;

    // Hadeler rejects exactly when all six distinct entries of the symmetric adjugate are positive.
    fraction::mul(work, b22, b33);
    work.submul(b23, b23);
    if (work.sgn() <= 0) return true;
    fraction::mul(work, b11, b33);
    work.submul(b13, b13);
    if (work.sgn() <= 0) return true;
    fraction::mul(work, b11, b22);
    work.submul(b12, b12);
    if (work.sgn() <= 0) return true;
    fraction::mul(work, b13, b23);
    work.submul(b12, b33);
    if (work.sgn() <= 0) return true;
    fraction::mul(work, b12, b23);
    work.submul(b13, b22);
    if (work.sgn() <= 0) return true;
    fraction::mul(work, b12, b13);
    work.submul(b11, b23);
    return work.sgn() <= 0;
}

/*
 * Exact strict-copositivity test for a symmetric rational matrix A:
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
    // Compute adj(A) from cofactors when A is singular and A^-1 is unavailable.
    matrix_frc adjugate(const matrix_frc& A) {
        const size_t n = A.rows();
        matrix_frc adj(n, n);
        matrix_frc minor(n - 1, n - 1);

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
                        minor(minor_row, minor_col) = A(row, col);
                        ++minor_col;
                    }
                    ++minor_row;
                }

                // Its signed determinant is cofactor C_ij.
                LU_Factorization lu(minor);
                fraction det_minor = lu.determinant();
                fraction cofactor = ((i + j) & 1) == 0 ? det_minor : -det_minor;

                // adj(A) is the transposed cofactor matrix; symmetry makes both
                // mirrored entries equal.
                adj(j, i) = cofactor;
                if (i != j) {
                    adj(i, j) = cofactor;  // Symmetry: same cofactor value
                }
            }
        }
        return adj;
    }

    bool all_entries_greater_zero(const matrix_frc& A) const noexcept {
        const size_t n = A.rows();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                // Only one triangle is needed because every input here is symmetric.
                if (A(i, j).sgn() <= 0) {
                    return false;
                }
            }
        }
        return true;
    }

    bool is_copositive_hadeler(const matrix_frc& A, bitset64 mask, size_t current_dim) {
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
        matrix_frc subMat(current_dim, current_dim);
        for (size_t row = 0; row < current_dim; ++row) {
            for (size_t column = 0; column < current_dim; ++column) {
                subMat(row, column) =
                    A(subset_indices[row], subset_indices[column]);
            }
        }

        LU_Factorization lu(subMat);
        fraction det = lu.determinant();

        if (det <= fraction::zero()) {
            /*
             * With all proper principal submatrices already passed, Hadeler's
             * rejecting pattern is det(A) <= 0 and adj(A) > 0 entrywise. In the
             * singular case, a positive adjugate gives a positive null direction,
             * so the quadratic form is not strictly positive there.
             */
            matrix_frc adj;
            if (lu.isSingular()) {
                // No inverse exists, so use the cofactor definition.
                adj = adjugate(subMat);
            } else {
                // For nonsingular A, adj(A) = det(A) * A^-1.
                adj = lu.inverse() * det;
            }

            if (all_entries_greater_zero(adj)) {
                return false;
            }
        }

        return true;
    }

public:
    bool is_strictly_copositive(const matrix_frc& A) {
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
};

inline bool is_strictly_copositive(const matrix_frc& A) {
    CopositivityChecker checker;
    return checker.is_strictly_copositive(A);
}

} // namespace linalg

#endif // RATIONAL_LINALG_COPOSITIVITY_FRACTION_HPP
