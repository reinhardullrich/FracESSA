#ifndef LINALG_COPOSITIVITY_INTEGER_HPP
#define LINALG_COPOSITIVITY_INTEGER_HPP

#include <fracessa/bitset64.hpp>
#include <fracessa/bitset_multiword.hpp>
#include <linalg/fraction_free_ldlt.hpp>
#include <linalg/matrix_integer.hpp>
#include <array>
#include <cassert>
#include <cstddef>
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

struct CopositivitySignScanMultiword {
    explicit CopositivitySignScanMultiword(size_t dimension)
    {
        negative_neighbors.reserve(dimension);
        for (size_t row = 0; row < dimension; ++row) negative_neighbors.emplace_back(dimension);
    }

    // As above, a failed diagonal leaves every other field intentionally incomplete.
    std::vector<bitset_multiword> negative_neighbors;
    integer all_ones_quadratic_value;
    bool all_diagonal_positive = true;
    bool has_negative_off_diagonal = false;
    bool all_negative_part_row_sums_positive = true;
};

// Large-dimension counterpart of scan_copositivity_signs(). It is separate so the one-word path keeps its fixed stack storage.
inline CopositivitySignScanMultiword scan_copositivity_signs_multiword(const matrix_int& A) {
    const size_t n = A.rows();
    CopositivitySignScanMultiword result(n);
    std::vector<integer> negative_part_row_sums(n);

    for (size_t i = 0; i < n; ++i) {
        const integer::const_reference diagonal = A(i, i);
        if (diagonal.sign() <= 0) {
            result.all_diagonal_positive = false;
            return result;
        }
        negative_part_row_sums[i] = diagonal;
        result.all_ones_quadratic_value += diagonal;
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const integer::const_reference entry = A(i, j);
            if (entry.sign() < 0) {
                result.negative_neighbors[i].set_bit_at_pos(j);
                result.negative_neighbors[j].set_bit_at_pos(i);
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
 * Principal submatrices are checked in increasing cardinality. When a matrix C is reached, every proper principal submatrix has
 * already passed. Hadeler (1983), Theorem 3, then decides C from its determinant, one retained solve, or its exact nullspace.
 */
class CopositivityChecker {
private:
    template<typename Index>
    bool decide_principal_submatrix_hadeler(const matrix_int& A, const Index* subset_indices, size_t current_dim) {
        matrix_int principal_submatrix(current_dim, current_dim);
        for (size_t row = 0; row < current_dim; ++row) {
            for (size_t column = 0; column < current_dim; ++column) {
                principal_submatrix(row, column) = A(subset_indices[row], subset_indices[column]);
            }
        }

        const bool nonsingular = factorization_.factorize_inplace(principal_submatrix) != 0;
        const int determinant_sign = factorization_.determinant().sign();
        if (determinant_sign > 0) return true;

        if (nonsingular) { // Here determinant_sign < 0.
            // Hadeler lets one positive right-hand side replace the complete adjugate: reject exactly when C y = -1 has y > 0.
            matrix_int solution(current_dim, 1);
            const integer minus_one(-1);
            for (size_t row = 0; row < current_dim; ++row) solution(row, 0) = minus_one;

            integer denominator;
            factorization_.solve_inplace(solution, denominator, principal_submatrix);
            assert(denominator.sign() > 0);

            for (size_t row = 0; row < current_dim; ++row) {
                if (solution(row, 0).sign() <= 0) return true;
            }
            return false;
        }

        // factorize_inplace() overwrote the matrix. Restore the original principal submatrix only for the singular nullspace test.
        for (size_t row = 0; row < current_dim; ++row) {
            for (size_t column = 0; column < current_dim; ++column) {
                principal_submatrix(row, column) = A(subset_indices[row], subset_indices[column]);
            }
        }

        matrix_int nullspace(current_dim, current_dim);
        const slong nullity = fmpz_mat_nullspace(nullspace.native_handle(), principal_submatrix.native_handle());
        assert(nullity > 0);
        if (nullity != 1) return true;

        // A one-dimensional nullspace rejects only when its basis is strictly one-sign; its orientation is arbitrary.
        const int basis_sign = nullspace(0, 0).sign();
        if (basis_sign == 0) return true;
        for (size_t row = 1; row < current_dim; ++row) {
            if (nullspace(row, 0).sign() != basis_sign) return true;
        }
        return false;
    }

    bool is_strictly_copositive_hadeler(const matrix_int& A, bitset64 mask, size_t current_dim) {
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
        return decide_principal_submatrix_hadeler(A, subset_indices, current_dim);
    }

    bool is_strictly_copositive_hadeler_multiword(const matrix_int& A, size_t current_dim) {
        switch (current_dim) {
            case 1: {
                const size_t i = subset_indices_[0];
                return is_strictly_copositive_1x1(A(i, i));
            }
            case 2: {
                const size_t i = subset_indices_[0];
                const size_t j = subset_indices_[1];
                return is_strictly_copositive_2x2(A(i, i), A(i, j), A(j, j));
            }
            case 3: {
                const size_t i = subset_indices_[0];
                const size_t j = subset_indices_[1];
                const size_t k = subset_indices_[2];
                return is_strictly_copositive_3x3(A(i, i), A(i, j), A(i, k), A(j, j), A(j, k), A(k, k));
            }
            default:
                return decide_principal_submatrix_hadeler(A, subset_indices_.data(), current_dim);
        }
    }

    bool enumerate_multiword_subsets(const matrix_int& A, size_t next_position, size_t needed, size_t subset_size) {
        if (needed == 0) return is_strictly_copositive_hadeler_multiword(A, subset_size);

        const size_t last_position = A.rows() - needed;
        // Recursion already knows each selected strategy, so retain its principal index instead of rebuilding indices from a mask.
        const size_t index_position = subset_size - needed;
        for (size_t position = next_position; position <= last_position; ++position) {
            subset_indices_[index_position] = position;
            if (!enumerate_multiword_subsets(A, position + 1, needed - 1, subset_size)) return false;
        }
        return true;
    }

    bool is_strictly_copositive_multiword(const matrix_int& A) {
        subset_indices_.resize(A.rows());
        for (size_t subset_size = 1; subset_size <= A.rows(); ++subset_size) {
            if (!enumerate_multiword_subsets(A, 0, subset_size, subset_size)) return false;
        }
        return true;
    }

public:
    explicit CopositivityChecker(size_t maximum_dimension)
        : factorization_(maximum_dimension)
    {
        if (maximum_dimension > bs64::kMaxBitsetDimension) subset_indices_.reserve(maximum_dimension);
    }

    bool is_strictly_copositive(const matrix_int& A) {
        const size_t n = A.rows();
        if (n > bs64::kMaxBitsetDimension) return is_strictly_copositive_multiword(A);

        // Cardinality order establishes Hadeler's proper-submatrix precondition. The last-subset sentinel also supports n=64 without
        // shifting a uint64_t by 64.
        for (size_t subset_size = 1; subset_size <= n; ++subset_size) {
            bitset64 subset = bs64::set_all_n_bits(subset_size);
            const bitset64 last_subset = subset << (n - subset_size);
            while (true) {
                if (!is_strictly_copositive_hadeler(A, subset, subset_size)) return false;
                if (subset == last_subset) break;
                subset = bs64::next_same_cardinality(subset);
            }
        }

        return true;
    }

    // Caller has already proved that every diagonal entry is positive and provides the negative-entry graph from the shared sign
    // scan. Entries between distinct graph components are nonnegative, so A is strictly copositive exactly when every component's
    // principal matrix is strictly copositive. Components are found once here; Hadeler does not split its principal subsets again.
    bool are_negative_components_strictly_copositive(
        const matrix_int& A, const std::array<bitset64, bs64::kMaxBitsetDimension>& negative_neighbors) {
        const size_t n = A.rows();
        const bitset64 all_vertices = bs64::set_all_n_bits(n);
        bitset64 remaining = all_vertices;

        while (remaining != 0) {
            bitset64 component = bs64::lowest_set_bit_as_bit(remaining);
            bitset64 frontier = component;
            while (frontier != 0) {
                const size_t vertex = bs64::find_pos_first_set_bit(frontier);
                frontier &= frontier - 1;
                const bitset64 discovered = negative_neighbors[vertex] & ~component;
                component |= discovered;
                frontier |= discovered;
            }

            // The common connected case needs no component matrix or allocation.
            if (component == all_vertices) return is_strictly_copositive(A);

            remaining &= ~component;
            if ((component & (component - 1)) == 0) continue;

            uint8_t indices[bs64::kMaxBitsetDimension];
            const size_t component_size = bs64::extract_set_indices(component, n, indices);
            matrix_int component_matrix(component_size, component_size);
            for (size_t row = 0; row < component_size; ++row) {
                for (size_t column = 0; column < component_size; ++column) {
                    component_matrix(row, column) = A(indices[row], indices[column]);
                }
            }
            if (!is_strictly_copositive(component_matrix)) return false;
        }

        return true;
    }

    bool are_negative_components_strictly_copositive(
        const matrix_int& A, const std::vector<bitset_multiword>& negative_neighbors) {
        const size_t n = A.rows();
        assert(n > bs64::kMaxBitsetDimension);
        assert(negative_neighbors.size() == n);

        bitset_multiword all_vertices(n);
        all_vertices.set_all();
        bitset_multiword remaining(all_vertices);
        bitset_multiword component(n);
        bitset_multiword frontier(n);
        bitset_multiword discovered(n);

        while (!remaining.empty()) {
            component.clear_all();
            frontier.clear_all();
            const size_t first_vertex = remaining.find_pos_first_set_bit();
            component.set_bit_at_pos(first_vertex);
            frontier.set_bit_at_pos(first_vertex);

            while (!frontier.empty()) {
                const size_t vertex = frontier.find_pos_first_set_bit();
                frontier.clear_bit_at_pos(vertex);
                discovered = negative_neighbors[vertex];
                discovered.subtract(component);
                component.union_with(discovered);
                frontier.union_with(discovered);
            }

            if (component == all_vertices) return is_strictly_copositive(A);

            remaining.subtract(component);
            const size_t component_size = component.count_set_bits();
            if (component_size == 1) continue;

            component.extract_set_indices(subset_indices_);
            matrix_int component_matrix(component_size, component_size);
            for (size_t row = 0; row < component_size; ++row) {
                for (size_t column = 0; column < component_size; ++column) {
                    component_matrix(row, column) = A(subset_indices_[row], subset_indices_[column]);
                }
            }
            if (!is_strictly_copositive(component_matrix)) return false;
        }

        return true;
    }

private:
    fraction_free_ldlt_factorization factorization_;
    std::vector<size_t> subset_indices_;
};

inline bool is_strictly_copositive(const matrix_int& A) {
    CopositivityChecker checker(A.rows());
    return checker.is_strictly_copositive(A);
}

} // namespace linalg

#endif // LINALG_COPOSITIVITY_INTEGER_HPP
