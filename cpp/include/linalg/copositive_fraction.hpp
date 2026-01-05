#ifndef RATIONAL_LINALG_COPOSITIVE_FRACTION_HPP
#define RATIONAL_LINALG_COPOSITIVE_FRACTION_HPP

#include <linalg/matrix_fraction.hpp>
#include <linalg/lu_factor_fraction.hpp>
#include <fracessa/bitset64.hpp>
#include <unordered_map>

namespace linalg {

inline matrix_frc adjugate(const matrix_frc& A) {
    const size_t n = A.rows();
    if (n == 0) return matrix_frc(0, 0);
    if (n == 1) {
        matrix_frc result(1, 1);
        result(0, 0) = fraction::one();
        return result;
    }
    
    // Hardcoded 2x2 adjugate for symmetric matrices (Bee matrices are symmetric)
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
    

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
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
            
            LU_Factorization lu(minor);
            fraction det_minor = lu.determinant();
            fraction val = ((i + j) & 1) == 0 ? det_minor : -det_minor;
            
            // Adjugate is the transpose of the cofactor matrix.
            // For symmetric A, C is symmetric, so transposing doesn't change it.
            // However, cofactor C(i,j) maps to adj(j,i).
            adj(j, i) = val;
            if (i != j) adj(i, j) = val;
        }
    }
    return adj;
}

struct bitset64_hash {
    std::size_t operator()(const bitset64& bs) const noexcept {
        return bs64::hash(bs);
    }
};

class CopositivityChecker {
private:
    std::unordered_map<bitset64, bool, bitset64_hash> memo;

    bool checkRecursive(const matrix_frc& A, const bitset64& mask) {
        auto it = memo.find(mask);
        if (it != memo.end()) return it->second;

        size_t current_dim = bs64::count_set_bits(mask);

        if (current_dim == 1) {
            size_t idx = bs64::find_pos_first_set_bit(mask);
            bool result = A(idx, idx) > fraction::zero();
            memo[mask] = result;
            return result;
        }

        // 3. Recursive Step: Check all proper principal submatrices of size (current_dim - 1)
        // A matrix is strictly copositive ONLY IF all its proper principal submatrices are.
        bool all_good = true;
        for (size_t i = bs64::find_pos_first_set_bit(mask); i < 64; i = bs64::find_pos_next_set_bit(mask, i)) {
            bitset64 sub_mask = bs64::clear_bit_at_pos(mask, i);
            if (!checkRecursive(A, sub_mask)) {
                all_good = false;
                break;
            }
        }
        if (!all_good) {
            memo[mask] = false;
            return false;
        }

        /**
         * 4. Hadeler / Motzkin-Straus / Kaplansky Criterion:
         * If all proper principal submatrices are strictly copositive,
         * A is strictly copositive UNLESS (det(A) <= 0 AND adj(A) has any non-positive entry).
         * 
         * More precisely: A is strictly copositive if and only if:
         *  a) All proper principal submatrices are strictly copositive.
         *  b) AND EITHER det(A) > 0 OR some entry of adj(A) is <= 0.
         */
        // Create submatrix in-place using reference logic as requested
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
            // Compute Adjugate to check if it's strictly positive
            matrix_frc adj;
            if (lu.isSingular()) {
                // Fall back to cofactor method for singular matrices
                // Since submatrix of symmetric A is symmetric, adjugate remains optimized.
                adj = adjugate(subMat);
            } else {
                // For det < 0, we can use the inverse formula: adj(A) = det(A) * A^(-1)
                adj = lu.inverse() * det;
            }

            if (adj.all_entries_greater_zero()) {
                memo[mask] = false;
                return false;
            }
        }

        memo[mask] = true;
        return true;
    }

public:
    bool isStrictlyCopositiveMemoized(const matrix_frc& A) {
        size_t n = A.rows();
        memo.clear();
        memo.reserve(bs64::two_to_the_power_of(n));
        bitset64 full_mask = bs64::set_all_n_bits(n);
        return checkRecursive(A, full_mask);
    }
};

inline bool isStrictlyCopositiveMemoized(const matrix_frc& A) {
    CopositivityChecker checker;
    return checker.isStrictlyCopositiveMemoized(A);
}

} // namespace linalg

#endif // RATIONAL_LINALG_COPOSITIVE_FRACTION_HPP
