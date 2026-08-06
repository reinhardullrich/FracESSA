#ifndef LINALG_COPOSITIVITY_INTEGER_HPP
#define LINALG_COPOSITIVITY_INTEGER_HPP

#include <fracessa/bitset64.hpp>
#include <linalg/matrix_integer.hpp>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace linalg {

/*
 * Facts obtained from one exact sign scan of a symmetric matrix. The negative-neighbor masks retain the graph needed by the later
 * connected-component reduction; the other fields answer the remaining cheap sign questions without another matrix pass.
 */
struct CopositivitySignScan {
    // Caller contract: test all_diagonal_positive before reading any other field. A false value means the scan returned before the
    // off-diagonal pass, so every other field is intentionally incomplete.
    std::array<bitset64, bs64::kMaxBitsetDimension> negative_neighbors{};
    integer all_ones_quadratic_value;
    bool all_diagonal_positive = true;
    bool has_negative_off_diagonal = false;
    bool all_negative_part_row_sums_positive = true;
};

/*
 * Scan the integer-scaled reduced B matrix once. The quadratic sum and every
 * sign decision remain exact even when entries exceed a machine word.
 */
inline CopositivitySignScan scan_copositivity_signs(const matrix_int& A) {
    CopositivitySignScan result;
    std::array<integer, 64> negative_part_row_sums;
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

    // Symmetry lets one triangular scan update both row sums and count every off-diagonal contribution to 1^T A 1 exactly twice.
    for (size_t i = 0; i < n; ++i) {
        const bitset64 bit_i = bs64::single_bit_at_pos(i);
        for (size_t j = i + 1; j < n; ++j) {
            const integer::const_reference entry = A(i, j);
            if (entry.sign() < 0) {
                const bitset64 bit_j = bs64::single_bit_at_pos(j);
                result.negative_neighbors[i] |= bit_j;
                result.negative_neighbors[j] |= bit_i;
                result.has_negative_off_diagonal = true;
                negative_part_row_sums[i] += entry;
                negative_part_row_sums[j] += entry;
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
 * FracESSA reaches this code only after the simpler exact low-dimensional and sign cases have failed. The remaining test adaptively
 * subdivides the nonnegative orthant into simplicial cones. Each cone is stored
 * only through B = V A V^T; replacing one generator v_i by v_i + v_j therefore
 * needs one exact row-and-column update instead of rebuilding the basis V.
 *
 * Algorithmic origin: Mathieu Dutour Sikiric's PairDecomposition and
 * TestStrictCopositivity in Polyhedral Common, src_copos/Copositivity.h,
 * commit d2252bc89d991fa6df9750ac9647e19b6a9aca02:
 * https://github.com/MathieuDutSik/polyhedral_common/blob/d2252bc89d991fa6df9750ac9647e19b6a9aca02/src_copos/Copositivity.h
 * This implementation was independently written for FracESSA. It uses FLINT
 * integers, retains no witness or basis, and updates the Gram matrix directly.
 */
class CopositivityChecker {
private:
    static bool inspect_cone(const matrix_int& matrix, size_t& split_i, size_t& split_j) {
        const size_t n = matrix.rows();
        split_i = n;
        split_j = n;

        for (size_t i = 0; i < n; ++i) {
            if (matrix(i, i).sign() <= 0) return false;
        }

        integer best_numerator;
        integer best_denominator;
        integer numerator;
        integer denominator;
        integer left;
        integer right;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                const integer::const_reference entry = matrix(i, j);
                if (entry.sign() >= 0) continue;

                numerator.set_product(entry, entry);
                denominator.set_product(matrix(i, i), matrix(j, j));
                const int local_comparison = numerator.compare(denominator);

                // Equality gives a zero direction; a larger numerator gives a negative direction in span{v_i, v_j}.
                if (local_comparison >= 0) return false;

                bool take_pair = split_i == n;
                if (!take_pair) {
                    left.set_product(numerator, best_denominator);
                    right.set_product(best_numerator, denominator);
                    take_pair = left.compare(right) > 0;
                }
                if (take_pair) {
                    split_i = i;
                    split_j = j;
                    best_numerator = numerator;
                    best_denominator = denominator;
                }
            }
        }

        return true;
    }

    static void replace_generator_with_sum(matrix_int& matrix, size_t replaced, size_t other) {
        integer new_diagonal(matrix(replaced, replaced));
        new_diagonal += matrix(replaced, other);
        new_diagonal += matrix(replaced, other);
        new_diagonal += matrix(other, other);

        const size_t n = matrix.rows();
        for (size_t k = 0; k < n; ++k) {
            if (k == replaced) continue;
            integer sum(matrix(replaced, k));
            sum += matrix(other, k);
            matrix(replaced, k) = sum;
            matrix(k, replaced) = sum;
        }
        matrix(replaced, replaced) = new_diagonal;
    }

public:
    static bool is_strictly_copositive(const matrix_int& A) {
        switch (A.rows()) {
            case 1:
                return is_strictly_copositive_1x1(A(0, 0));
            case 2:
                return is_strictly_copositive_2x2(A(0, 0), A(0, 1), A(1, 1));
            case 3:
                return is_strictly_copositive_3x3(A(0, 0), A(0, 1), A(0, 2), A(1, 1), A(1, 2), A(2, 2));
            default:
                break;
        }

        std::vector<matrix_int> pending;
        pending.reserve(64);
        pending.emplace_back(A);

        while (!pending.empty()) {
            matrix_int current(std::move(pending.back()));
            pending.pop_back();

            size_t split_i;
            size_t split_j;
            if (!inspect_cone(current, split_i, split_j)) return false;
            if (split_i == current.rows()) continue;

            // Retain one child in current and copy the unsplit Gram matrix only once for its sibling.
            matrix_int sibling(current);
            replace_generator_with_sum(current, split_j, split_i);
            replace_generator_with_sum(sibling, split_i, split_j);
            pending.emplace_back(std::move(current));
            pending.emplace_back(std::move(sibling));
        }

        return true;
    }
};

} // namespace linalg

#endif // LINALG_COPOSITIVITY_INTEGER_HPP
