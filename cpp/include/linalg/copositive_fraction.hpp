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
 * Strict copositivity checker for exact rational symmetric matrices.
 *
 * Context in FracESSA:
 * - Called only after cheaper ESS filters fail (pure ESS / PD checks).
 * - Works on Bee-type reduced matrices where exact sign decisions matter.
 *
 * Method sketch:
 * 1) Enumerate principal submatrices by support subsets.
 * 2) For each subset, evaluate the Hadeler/Motzkin-Straus/Kaplansky criterion.
 * 3) Abort early on the first violating subset.
 */
class CopositivityCheckerV3 {
private:
    std::vector<std::vector<bitset64>> supports_;

    // Exact adjugate. Needed when LU inverse is unavailable (singular case),
    // but the criterion still requires checking adj(A) entry signs.
    matrix_frc adjugate(const matrix_frc& A) {
        const size_t n = A.rows();
        if (n == 0) return matrix_frc(0, 0);
        if (n == 1) {
            matrix_frc result(1, 1);
            result(0, 0) = fraction::one();
            return result;
        }
        // Hardcoded 2x2 adjugate for symmetric matrices (Bee matrices are symmetric, bee_vee matrices as well!)
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

        // Optimization: For symmetric A, we only need to compute cofactors for upper triangle
        // since C(i,j) = C(j,i) for symmetric matrices (det(minor(i,j)) = det(minor(j,i)))
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                // Extract minor(i,j) - matrix with row i and column j removed
                size_t minor_row = 0;
                for (size_t row = 0; row < n; ++row) {
                    if (row == i) continue;
                    size_t minor_col = 0;
                    for (size_t col = 0; col < n; ++col) {
                        if (col == j) continue;
                        // Use symmetry: A(row, col) = A(col, row) when extracting
                        minor(minor_row, minor_col) = A(row, col);
                        ++minor_col;
                    }
                    ++minor_row;
                }

                // Compute determinant once - for symmetric A, det(minor(i,j)) = det(minor(j,i))
                LU_Factorization lu(minor);
                fraction det_minor = lu.determinant();
                fraction cofactor = ((i + j) & 1) == 0 ? det_minor : -det_minor;

                // Adjugate is the transpose of the cofactor matrix
                // For symmetric A: adj(j,i) = C(i,j) and adj(i,j) = C(j,i) = C(i,j)
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
                // Use fraction's sgn() method which internally uses fmpq_sgn
                // This avoids creating temporary fraction objects for comparison
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

        // Build principal submatrix A[mask,mask] in compact coordinates.
        matrix_frc subMat(current_dim, current_dim);
        {
            size_t row = 0;
            for (size_t i = bs64::find_pos_first_set_bit(mask); i < 64; i = bs64::find_pos_next_set_bit(mask, i)) {
                size_t col = 0;
                for (size_t j = bs64::find_pos_first_set_bit(mask); j < 64; j = bs64::find_pos_next_set_bit(mask, j)) {
                    // Using direct assignment. Since LU_Factorization copies anyway,
                    // we avoid double initialization by using default constructor (no init).
                    subMat(row, col) = A(i, j);
                    ++col;
                }
                ++row;
            }
        }

        LU_Factorization lu(subMat);
        fraction det = lu.determinant();

        if (det <= fraction::zero()) {
            /**
             * 4. Hadeler / Motzkin-Straus / Kaplansky Criterion:
             * If all proper principal submatrices are strictly copositive,
             * A is strictly copositive UNLESS (det(A) <= 0 AND adj(A) has any non-positive entry).
             *
             * More precisely: A is strictly copositive if and only if:
             *  a) All proper principal submatrices are strictly copositive.
             *  b) AND EITHER det(A) > 0 OR some entry of adj(A) is <= 0.
             */
            // Compute adj(A) for the second part of the criterion.
            matrix_frc adj;
            if (lu.isSingular()) {
                // Fall back to cofactor method for singular matrices
                // Since submatrix of symmetric A is symmetric, adjugate remains optimized.
                adj = adjugate(subMat);
            } else {
                // For det < 0, we can use the inverse formula: adj(A) = det(A) * A^(-1)
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

        // Group all non-empty supports by cardinality once.
        supports_.resize(n);
        for (size_t i = 0; i < n; ++i) {
            supports_[i].reserve(binomial_coefficient(n, i + 1));
        }
        for (uint64_t bits = 1ULL; bits < bs64::two_to_the_power_of(n); ++bits) {
            bitset64 support = bits;
            supports_[bs64::count_set_bits(support) - 1].push_back(support);
        }

        // Check from small subsets to full matrix. Any failure is conclusive.
        for (size_t subset_size = 1; subset_size <= n; ++subset_size) {
            const auto& subsets_of_size = supports_[subset_size - 1];
            for (bitset64 subset : subsets_of_size) {
                if (!is_copositive_hadeler(A, subset, subset_size)) {
                    return false;  // Early exit - found non-copositive submatrix
                }
            }
        }

        // Passed all subset checks: matrix is strictly copositive.
        return true;
    }
};

inline bool is_strictly_copositive(const matrix_frc& A) {
    CopositivityCheckerV3 checker;
    return checker.is_strictly_copositive(A);
}

} // namespace linalg

#endif // RATIONAL_LINALG_COPOSITIVITY_FRACTION_HPP
