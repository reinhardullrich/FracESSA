#ifndef RATIONAL_LINALG_COPOSITIVITY_FRACTION_HPP
#define RATIONAL_LINALG_COPOSITIVITY_FRACTION_HPP

#include <linalg/matrix_fraction.hpp>
#include <linalg/lu_factor_fraction.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/supports.hpp>
#include <cstddef>
#include <vector>

namespace linalg {

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
class CopositivityCheckerV3 {
private:
    std::vector<std::vector<bitset64>> supports_;

    // Compute adj(A) from cofactors when A is singular and A^-1 is unavailable.
    matrix_frc adjugate(const matrix_frc& A) {
        const size_t n = A.rows();
        if (n == 0) return matrix_frc(0, 0);
        if (n == 1) {
            matrix_frc result(1, 1);
            result(0, 0) = fraction::one();
            return result;
        }
        // The common 2-by-2 case is clearer and cheaper without general minors.
        if (n == 2) {
            matrix_frc result(2, 2);
            result(0, 0) = A(1, 1);
            result(1, 1) = A(0, 0);
            result(0, 1) = -A(0, 1); // Symmetric, so A(0,1) == A(1,0)
            result(1, 0) = result(0, 1);
            return result;
        }

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

        if (current_dim == 1) {
            size_t idx = bs64::find_pos_first_set_bit(mask);
            return A(idx, idx) > fraction::zero();
        }

        // Copy the principal submatrix A[mask,mask] into compact row/column order.
        matrix_frc subMat(current_dim, current_dim);
        {
            size_t row = 0;
            for (size_t i = bs64::find_pos_first_set_bit(mask); i < 64; i = bs64::find_pos_next_set_bit(mask, i)) {
                size_t col = 0;
                for (size_t j = bs64::find_pos_first_set_bit(mask); j < 64; j = bs64::find_pos_next_set_bit(mask, j)) {
                    subMat(row, col) = A(i, j);
                    ++col;
                }
                ++row;
            }
        }

        LU_Factorization lu(subMat);
        fraction det = lu.determinant();

        if (det <= fraction::zero()) {
            /*
             * With all proper principal submatrices already passed, Hadeler's
             * rejecting pattern is det(A) < 0 and adj(A) > 0 entrywise. The code
             * enters for det(A) == 0 as well because a singular matrix has no
             * inverse; under the theorem's precondition its adjugate cannot
             * match the rejecting pattern.
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

        // Each nonempty mask selects one principal submatrix. Grouping by size
        // establishes the proper-submatrix precondition used by the criterion.
        supports_.resize(n);
        for (size_t i = 0; i < n; ++i) {
            supports_[i].reserve(binomial_coefficient(n, i + 1));
        }
        for (uint64_t bits = 1ULL; bits < bs64::two_to_the_power_of(n); ++bits) {
            bitset64 support = bits;
            supports_[bs64::count_set_bits(support) - 1].push_back(support);
        }

        // Stop at the first principal submatrix that is not strictly copositive.
        for (size_t subset_size = 1; subset_size <= n; ++subset_size) {
            const auto& subsets_of_size = supports_[subset_size - 1];
            for (bitset64 subset : subsets_of_size) {
                if (!is_copositive_hadeler(A, subset, subset_size)) {
                    return false;
                }
            }
        }

        // In particular, the final full-size subset passed.
        return true;
    }
};

inline bool is_strictly_copositive(const matrix_frc& A) {
    CopositivityCheckerV3 checker;
    return checker.is_strictly_copositive(A);
}

} // namespace linalg

#endif // RATIONAL_LINALG_COPOSITIVITY_FRACTION_HPP
