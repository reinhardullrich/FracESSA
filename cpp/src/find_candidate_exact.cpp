#include <fracessa/find_candidate_exact.hpp>

#include <array>
#include <cstdint>
#include <utility>

namespace candidate_search {
namespace {

using permutation_array = std::array<uint8_t, bs64::kMaxBitsetDimension>;

struct ldlt_result {
    bool nonsingular;
    bool negative_definite;
};

/*
 * The exact LDL^T factorization below stores only the lower triangle:
 *
 *   - entries strictly below a completed pivot block belong to the unit-lower triangular matrix L;
 *   - a 1x1 block of D is stored on the diagonal;
 *   - a 2x2 block of D is stored in its two diagonal positions and the one lower off-diagonal position.
 *
 * The upper triangle is not read. Accessing the active symmetric Schur complement through this helper therefore keeps every update in one place.
 */
fraction& lower_entry(linalg::matrix_frc& matrix, size_t row, size_t column)
{
    return row >= column ? matrix(row, column) : matrix(column, row);
}

/*
 * Apply one symmetric permutation to the partially factored matrix.
 *
 * A row-only swap would destroy symmetry and would not preserve inertia. The same two indices must therefore be exchanged as both rows and columns.
 * The already completed columns contain L rather than the Schur complement, so their two row entries are swapped separately. The active lower
 * triangle is then permuted as a symmetric matrix, and the right-hand side receives the matching row permutation.
 */
void swap_symmetric_indices(linalg::matrix_frc& system, size_t dimension, size_t completed_columns, size_t first, size_t second)
{
    if (first == second) return;

    for (size_t column = 0; column < completed_columns; ++column) {
        std::swap(system(first, column), system(second, column));
    }

    std::swap(system(first, dimension), system(second, dimension));
    std::swap(system(first, first), system(second, second));

    // The entry joining first and second is unchanged by exchanging both of its indices. Every other active entry touching either index is swapped.
    for (size_t index = completed_columns; index < dimension; ++index) {
        if (index == first || index == second) continue;
        std::swap(lower_entry(system, first, index), lower_entry(system, second, index));
    }
}

/*
 * Factor and solve one exact symmetric reduced system in place.
 *
 * On entry, the lower n-by-n triangle is the reduced Hessian H and column n is the right-hand side r. On successful return, that column contains
 * the solution in the current permuted order, while permutation[position] gives the original reduced coordinate represented by that position.
 *
 * Exact arithmetic removes the numerical reason for Bunch-Kaufman threshold choices, but symmetric indefinite matrices still need pivoting for
 * algebraic correctness. At each step:
 *
 *   1. Any nonzero diagonal entry is a valid 1x1 pivot.
 *   2. If every remaining diagonal is zero, any nonzero off-diagonal entry forms the nonsingular block [[0,a],[a,0]].
 *   3. If neither exists, the remaining Schur complement is the zero matrix, so H is singular.
 *
 * These symmetric permutations and block eliminations are congruences. By Sylvester's law of inertia, H is negative definite exactly when every
 * 1x1 D pivot is negative and no 2x2 block occurs. A zero-diagonal 2x2 block always has one positive and one negative eigenvalue.
 */
ldlt_result factor_and_solve_reduced_system(linalg::matrix_frc& system, size_t dimension, permutation_array& permutation)
{
    // A value of 1 or 2 marks the first position of a D block; zero marks the second position of a 2x2 block.
    std::array<uint8_t, bs64::kMaxBitsetDimension> block_sizes{};
    for (size_t i = 0; i < dimension; ++i) {
        permutation[i] = static_cast<uint8_t>(i);
    }

    bool negative_definite = true;
    size_t pivot_position = 0;

    while (pivot_position < dimension) {
        // Prefer a 1x1 pivot. Magnitude comparisons are unnecessary because every operation and every zero test is exact.
        size_t diagonal_pivot = dimension;
        for (size_t i = pivot_position; i < dimension; ++i) {
            if (!system(i, i).is_zero()) {
                diagonal_pivot = i;
                break;
            }
        }

        if (diagonal_pivot != dimension) {
            swap_symmetric_indices(system, dimension, pivot_position, pivot_position, diagonal_pivot);
            std::swap(permutation[pivot_position], permutation[diagonal_pivot]);

            block_sizes[pivot_position] = 1;
            const fraction pivot = system(pivot_position, pivot_position);
            if (pivot.sgn() >= 0) negative_definite = false;

            // L(i,k)=H(i,k)/D(k,k). Compute every multiplier before the trailing Schur complement overwrites its source column.
            for (size_t row = pivot_position + 1; row < dimension; ++row) {
                system(row, pivot_position).div_inplace(pivot);
            }

            // H_next = H_trailing - L(:,k)*D(k,k)*L(:,k)^T.
            fraction scaled_column_entry;
            for (size_t row = pivot_position + 1; row < dimension; ++row) {
                for (size_t column = pivot_position + 1; column <= row; ++column) {
                    fraction::mul(scaled_column_entry, pivot, system(column, pivot_position));
                    system(row, column).submul(system(row, pivot_position), scaled_column_entry);
                }
            }

            ++pivot_position;
            continue;
        }

        // Every remaining diagonal is zero. A nonzero off-diagonal entry is therefore the only possible nonsingular symmetric pivot.
        size_t first_pivot = dimension;
        size_t second_pivot = dimension;
        for (size_t row = pivot_position + 1; row < dimension && first_pivot == dimension; ++row) {
            for (size_t column = pivot_position; column < row; ++column) {
                if (!system(row, column).is_zero()) {
                    first_pivot = column;
                    second_pivot = row;
                    break;
                }
            }
        }

        if (first_pivot == dimension) {
            return {false, false};
        }

        swap_symmetric_indices(system, dimension, pivot_position, pivot_position, first_pivot);
        std::swap(permutation[pivot_position], permutation[first_pivot]);

        swap_symmetric_indices(system, dimension, pivot_position, pivot_position + 1, second_pivot);
        std::swap(permutation[pivot_position + 1], permutation[second_pivot]);

        block_sizes[pivot_position] = 2;
        block_sizes[pivot_position + 1] = 0;
        negative_definite = false;

        // The selected D block is [[0,a],[a,0]], whose inverse is [[0,1/a],[1/a,0]]. Multiplication by that inverse exchanges the two entries
        // in every active row and divides both by a.
        const fraction off_diagonal_pivot = system(pivot_position + 1, pivot_position);
        for (size_t row = pivot_position + 2; row < dimension; ++row) {
            const fraction first_entry = system(row, pivot_position);
            const fraction second_entry = system(row, pivot_position + 1);
            fraction::div(system(row, pivot_position), second_entry, off_diagonal_pivot);
            fraction::div(system(row, pivot_position + 1), first_entry, off_diagonal_pivot);
        }

        // Subtract L_block*D_block*L_block^T from the trailing matrix. Because D has only the two off-diagonal entries a, each update has exactly
        // the two cross-products written below.
        fraction scaled_column_entry;
        for (size_t row = pivot_position + 2; row < dimension; ++row) {
            for (size_t column = pivot_position + 2; column <= row; ++column) {
                fraction::mul(scaled_column_entry, off_diagonal_pivot, system(column, pivot_position + 1));
                system(row, column).submul(system(row, pivot_position), scaled_column_entry);

                fraction::mul(scaled_column_entry, off_diagonal_pivot, system(column, pivot_position));
                system(row, column).submul(system(row, pivot_position + 1), scaled_column_entry);
            }
        }

        pivot_position += 2;
    }

    /*
     * Solve L*z=P^T*r by forward substitution. The lower off-diagonal entry inside a 2x2 pivot belongs to D, not L, and must be skipped. The
     * right-hand side column is overwritten with z because its earlier entries are never needed again in their original form.
     */
    for (size_t row = 0; row < dimension; ++row) {
        for (size_t column = 0; column < row; ++column) {
            if (block_sizes[column] == 2 && row == column + 1) continue;
            system(row, dimension).submul(system(row, column), system(column, dimension));
        }
    }

    // Solve D*w=z one pivot block at a time, again in the right-hand side column. Every pivot is known to be nonsingular from the factorization.
    for (size_t block = 0; block < dimension; ) {
        if (block_sizes[block] == 1) {
            system(block, dimension).div_inplace(system(block, block));
            ++block;
            continue;
        }

        const fraction first_rhs = system(block, dimension);
        const fraction second_rhs = system(block + 1, dimension);
        const fraction& off_diagonal_pivot = system(block + 1, block);
        fraction::div(system(block, dimension), second_rhs, off_diagonal_pivot);
        fraction::div(system(block + 1, dimension), first_rhs, off_diagonal_pivot);
        block += 2;
    }

    // Solve L^T*y=w backwards. As in the forward solve, skip the one stored D entry inside every 2x2 block. The result remains in permuted coordinates.
    for (size_t row = dimension; row-- > 0; ) {
        for (size_t column = row + 1; column < dimension; ++column) {
            if (block_sizes[row] == 2 && column == row + 1) continue;
            system(row, dimension).submul(system(column, row), system(column, dimension));
        }
    }

    return {true, negative_definite};
}

} // namespace

/*
 * Exact candidate construction for a proposed support S.
 *
 * A mixed strategy x with support S is a symmetric Nash equilibrium when
 *
 *   A_S x = u*1,          every used strategy earns the same payoff u,
 *   1^T x = 1,            probabilities sum to one,
 *   x_i > 0 for i in S,   S is the actual support,
 *   (A x)_i <= u outside S.
 */
find_candidate_exact::find_candidate_exact(const linalg::matrix_frc& game_matrix) noexcept
    : game_frc_(game_matrix)
    , dimension_(game_matrix.rows())
{
}

bool find_candidate_exact::find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
{
    reduced_hessian_is_negative_definite_ = false;

    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    const size_t support_count = bs64::extract_set_indices(support, dimension_, support_indices);
    const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
    const size_t non_support_count = bs64::extract_set_indices(complement, dimension_, non_support_indices);

    if (support_solution_.rows() != support_size) {
        support_solution_ = linalg::matrix_frc(support_size, 1);
    }

    /*
     * Eliminate the border before factorization.
     *
     * Let m be the first strategy in S and let Z have columns e_i-e_m for all other i in S. Every normalized support vector has the unique form
     *
     *     x = e_m + Z*y.
     *
     * Multiplying A_S*x=u*1 by Z^T eliminates the unknown payoff u and gives
     *
     *     H*y = r,       H = Z^T*A_S*Z,       r = -Z^T*A_S*e_m.
     *
     * In entries, for i,j in S without m,
     *
     *     H(i,j) = A(i,j)-A(i,m)-A(m,j)+A(m,m),
     *     r(i)   = A(m,m)-A(i,m).
     *
     * H has size (|S|-1)-by-(|S|-1), is symmetric, and is nonsingular exactly when the original bordered candidate matrix is nonsingular. Thus the
     * reduction loses neither solutions nor the singularity decision.
     */
    const size_t reference_index = support_indices[0];
    const size_t reduced_dimension = support_size - 1;

    if (reduced_dimension == 0) {
        // The tangent space of a pure support has dimension zero. Its reduced Hessian is therefore vacuously negative definite, normalization fixes
        // the sole probability at one, and no factorization is necessary.
        support_solution_(0, 0) = fraction::one();
        reduced_hessian_is_negative_definite_ = true;
    } else {
        if (reduced_system_.rows() != reduced_dimension) {
            // The final column holds r and is reused for all three triangular solves. The lower square triangle holds H and then its LDL^T.
            reduced_system_ = linalg::matrix_frc(reduced_dimension, reduced_dimension + 1);
        }

        const fraction& reference_diagonal = game_frc_(reference_index, reference_index);
        for (size_t row = 0; row < reduced_dimension; ++row) {
            const size_t i = support_indices[row + 1];

            reduced_system_(row, reduced_dimension) = reference_diagonal;
            reduced_system_(row, reduced_dimension) -= game_frc_(i, reference_index);

            // Only the lower triangle is needed by the symmetric factorization.
            for (size_t column = 0; column <= row; ++column) {
                const size_t j = support_indices[column + 1];
                fraction& value = reduced_system_(row, column);
                value = game_frc_(i, j);
                value -= game_frc_(i, reference_index);
                value -= game_frc_(reference_index, j);
                value += reference_diagonal;
            }
        }

        permutation_array permutation{};
        const ldlt_result factorization = factor_and_solve_reduced_system(reduced_system_, reduced_dimension, permutation);
        if (!factorization.nonsingular) return false;
        reduced_hessian_is_negative_definite_ = factorization.negative_definite;

        /*
         * The solved y values are in the symmetric permutation order chosen by LDL^T. Put them back into the original support order, then reconstruct
         * the eliminated reference probability from normalization:
         *
         *     x_m = 1 - sum(y).
         */
        fraction reference_probability = fraction::one();
        for (size_t position = 0; position < reduced_dimension; ++position) {
            const size_t original_reduced_position = permutation[position];
            const fraction& probability = reduced_system_(position, reduced_dimension);
            if (probability.sgn() <= 0) return false;

            support_solution_(original_reduced_position + 1, 0) = probability;
            reference_probability -= probability;
        }
        if (reference_probability.sgn() <= 0) return false;
        support_solution_(0, 0) = std::move(reference_probability);
    }

    // The reduction eliminated u. Recover it from the reference row of A_S*x; the reduced equations prove that every other support row has the same sum.
    fraction payoff = fraction::zero();
    for (size_t position = 0; position < support_count; ++position) {
        payoff.addmul(game_frc_(reference_index, support_indices[position]), support_solution_(position, 0));
    }

    result.extended_support = support;

    // Exact outside-support inequalities finish the Nash equilibrium check and collect every unused pure strategy tying the equilibrium payoff.
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        fraction rowsum = fraction::zero();
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            rowsum.addmul(game_frc_(i, support_indices[j_pos]), support_solution_(j_pos, 0));
        }

        if (rowsum > payoff) return false;
        if (rowsum == payoff) {
            result.extended_support = bs64::set_bit_at_pos(result.extended_support, i);
        }
    }

    result.payoff = payoff;
    result.payoff_dbl = payoff.to_dbl();
    result.extended_support_size = bs64::count_set_bits(result.extended_support);

    // Stability uses only the support, extended support, payoff, and the stored reduced-Hessian sign. Materialize the dense n-vector only for requested
    // candidate output or logging.
    if (materialize_vector) {
        if (result.vector.rows() != dimension_ || result.vector.cols() != 1) {
            result.vector = linalg::matrix_frc(dimension_, 1);
        }
        for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
            result.vector(non_support_indices[i_pos], 0) = fraction::zero();
        }
        for (size_t i_pos = 0; i_pos < support_count; ++i_pos) {
            result.vector(support_indices[i_pos], 0) = support_solution_(i_pos, 0);
        }
    }
    return true;
}

} // namespace candidate_search
