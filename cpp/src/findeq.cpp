#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/linear_solver.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

/*
 * Candidate construction for a proposed support S.
 *
 * A mixed strategy x with support S is a symmetric Nash equilibrium when
 *
 *   A_S x = u*1,          every used strategy earns the same payoff u,
 *   1^T x = 1,            probabilities sum to one,
 *   x_i > 0 for i in S,   S is the actual support,
 *   (A x)_i <= u outside S.
 *
 * The bordered linear system supplies x and u. Exact equality in the final
 * inequality identifies an additional best reply; these tied strategies form
 * the extended support J(x). The double routine may cheaply reject a support,
 * but only the rational routine can create a candidate.
 */

namespace {

// The danger estimate below is written specifically for IEEE-754 binary64.
static_assert(std::numeric_limits<double>::is_iec559);
static_assert(std::numeric_limits<double>::radix == 2);
static_assert(std::numeric_limits<double>::digits == 53);

// Extract S and its complement once. The nested matrix loops can then walk
// compact index arrays instead of testing one bit for every matrix entry.
inline void extract_support_partition(
    bitset64 support,
    size_t dimension,
    uint8_t (&support_indices)[bs64::kMaxBitsetDimension],
    size_t& support_count,
    uint8_t (&non_support_indices)[bs64::kMaxBitsetDimension],
    size_t& non_support_count) noexcept
{
    support_count = bs64::extract_set_indices(support, dimension, support_indices);
    non_support_count = 0;

    size_t support_pos = 0;
    for (size_t i = 0; i < dimension; ++i) {
        if (support_pos < support_count && i == static_cast<size_t>(support_indices[support_pos])) {
            ++support_pos;
        } else {
            non_support_indices[non_support_count++] = static_cast<uint8_t>(i);
        }
    }
}

} // namespace

bool fracessa::find_candidate_dbl(const bitset64& support, size_t support_size)
{
    /*
     * Unsafe rejection filter.
     *
     * `false` means the floating result looks clearly incompatible with a
     * candidate, so the support is skipped. `true` means "ask the exact solver";
     * it is returned for plausible candidates and for every numerically
     * suspicious case. This heuristic can still reject a valid candidate; exact
     * mode therefore bypasses it completely.
     */
    const auto& game_matrix = matrix_server_.get_game_matrix_dbl();
    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    size_t support_count = 0;
    size_t non_support_count = 0;
    extract_support_partition(
        support, dimension_, support_indices, support_count,
        non_support_indices, non_support_count);

    // Build [A'_S -1 | 0; 1^T 0 | 1] directly in reusable scratch. The exact
    // affine transformation preserves the game; conversion to A' remains heuristic.
    const size_t system_dimension = support_size + 1;
    auto& Ab = matrix_server_.get_linear_system_dbl(support_size);
    for (size_t row = 0; row < support_size; ++row) {
        const size_t i = support_indices[row];
        for (size_t column = 0; column < support_size; ++column) {
            Ab(row, column) = game_matrix(i, support_indices[column]);
        }
        Ab(row, support_size) = -1.0;
        Ab(row, system_dimension) = 0.0;
    }
    for (size_t column = 0; column < support_size; ++column) {
        Ab(support_size, column) = 1.0;
    }
    Ab(support_size, support_size) = 0.0;
    Ab(support_size, system_dimension) = 1.0;

    constexpr double kPivotTolerance = 1e-12;
    double minimum_pivot = std::numeric_limits<double>::infinity();
    for (size_t column = 0; column + 1 < system_dimension; ++column) {
        size_t max_row = column;
        double max_value = std::abs(Ab(column, column));
        for (size_t row = column + 1; row < system_dimension; ++row) {
            const double value = std::abs(Ab(row, column));
            if (value > max_value) {
                max_value = value;
                max_row = row;
            }
        }

        // A tiny or invalid pivot makes the floating proposal unreliable. It is
        // not a reason to reject the support, so hand it to exact arithmetic.
        if (!std::isfinite(max_value) || max_value < kPivotTolerance) return true;
        minimum_pivot = std::min(minimum_pivot, max_value);
        if (max_row != column) Ab.swap_rows(column, max_row);

        const double pivot = Ab(column, column);
        for (size_t row = column + 1; row < system_dimension; ++row) {
            const double factor = Ab(row, column) / pivot;
            for (size_t j = column + 1; j <= system_dimension; ++j) {
                Ab(row, j) -= factor * Ab(column, j);
            }
            Ab(row, column) = 0.0;
        }
    }

    const double last_pivot = std::abs(Ab(system_dimension - 1, system_dimension - 1));
    if (!std::isfinite(last_pivot) || last_pivot < kPivotTolerance) return true;
    minimum_pivot = std::min(minimum_pivot, last_pivot);

    double solution[bs64::kMaxBitsetDimension];
    double solution_scale = 1.0;
    for (size_t i = system_dimension; i-- > 0;) {
        double sum = Ab(i, system_dimension);
        for (size_t j = i + 1; j < system_dimension; ++j) {
            sum -= Ab(i, j) * solution[j];
        }
        const double value = sum / Ab(i, i);
        if (!std::isfinite(value)) return true;
        solution[i] = value;
        solution_scale = std::max(solution_scale, std::abs(value));
    }

    constexpr double kProbabilityTolerance = -1e-10;
    constexpr double kGuardMultiplier = 64.0;
    constexpr double kBinary64Epsilon = std::numeric_limits<double>::epsilon();
    const double base_risk =
        kGuardMultiplier * static_cast<double>(system_dimension) * kBinary64Epsilon;
    /*
     * Cheap danger veto. `margin` is the size of an apparent violation and
     * `minimum_pivot` is a rough measure of how strongly elimination could have
     * amplified rounding. If the estimated rounding risk is not clearly below
     * their product, the violation is suspicious and the exact solver decides.
     * This comparison is deliberately a heuristic, not a certified error bound.
     */
    const auto is_suspicious = [&](double margin, double decision_scale) noexcept {
        const double risk = base_risk * decision_scale;
        const double separation = minimum_pivot * margin;
        return !std::isfinite(risk) || !std::isfinite(separation) ||
               separation <= 0.0 || risk >= separation;
    };

    // Only a clearly negative support probability can justify an unsafe rejection.
    for (size_t i = 0; i < support_size; ++i) {
        if (solution[i] < kProbabilityTolerance) {
            return is_suspicious(-solution[i], solution_scale);
        }
    }

    const double payoff = solution[support_size];
    const double proposal_threshold = 1e-4 * static_cast<double>(dimension_);
    const double outside_scale =
        static_cast<double>(system_dimension) * solution_scale;
    if (!std::isfinite(outside_scale)) return true;
    // A strategy outside S invalidates the candidate if it earns more than u.
    // Small apparent gains go exact; only a separated gain reaches the veto.
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        double rowsum = 0.0;
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            rowsum += game_matrix(i, support_indices[j_pos]) * solution[j_pos];
        }
        const double gain = rowsum - payoff;
        if (!std::isfinite(gain)) return true;
        if (gain > proposal_threshold) {
            return is_suspicious(gain, outside_scale);
        }
    }

    return true;
}

bool fracessa::find_candidate_frc(const bitset64& support, size_t support_size)
{
    // Authoritative candidate test: every solve and comparison below is exact.
    const auto& game_matrix = matrix_server_.get_game_matrix_frc();
    auto& Ab = matrix_server_.get_linear_system_frc(support, support_size);
    
    linalg::matrix_frc solution;
    bool solved = linalg::solve_linear_frc(Ab, solution);
    
    if (!solved) return false;

    // Reuse the same branch-free partition strategy as the double path.
    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    size_t support_count = 0;
    size_t non_support_count = 0;
    extract_support_partition(support, dimension_, support_indices, support_count, non_support_indices, non_support_count);

    const fraction payoff = solution(support_size, 0);
    candidate_.extended_support = support;

    if (candidate_.vector.rows() != dimension_ || candidate_.vector.cols() != 1) {
        candidate_.vector = linalg::matrix_frc(dimension_, 1);
    }
    // Store the compact solution as a full game strategy for output and stability work.
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        candidate_.vector(non_support_indices[i_pos], 0) = fraction::zero();
    }
    for (size_t i_pos = 0; i_pos < support_count; ++i_pos) {
        candidate_.vector(support_indices[i_pos], 0) = solution(i_pos, 0);
    }
    
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        fraction rowsum = fraction::zero();
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            const size_t j = support_indices[j_pos];
            rowsum.addmul(game_matrix(i, j), solution(j_pos, 0));
        }
        
        // A profitable pure deviation proves that this is not a Nash equilibrium.
        if (rowsum > payoff) return false;
        
        // Equality means this unused strategy is also a best reply and belongs to J(x).
        if (rowsum == payoff)
            candidate_.extended_support = bs64::set_bit_at_pos(candidate_.extended_support, i);
    }

    candidate_.payoff = payoff;
    candidate_.payoff_dbl = payoff.to_dbl();        
    candidate_.extended_support_size = bs64::count_set_bits(candidate_.extended_support);
    
    return true;
}
