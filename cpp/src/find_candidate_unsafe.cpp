#include <fracessa/find_candidate_unsafe.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace candidate_search {
namespace {

// The danger estimate below is written specifically for IEEE-754 binary64.
static_assert(std::numeric_limits<double>::is_iec559);
static_assert(std::numeric_limits<double>::radix == 2);
static_assert(std::numeric_limits<double>::digits == 53);

} // namespace

find_candidate_unsafe::find_candidate_unsafe(const linalg::matrix_frc& game_matrix) noexcept
    : game_frc_(game_matrix)
    , dimension_(game_matrix.rows())
{
}

bool find_candidate_unsafe::normalize_game_matrix()
{
    const fraction& translation = game_frc_(0, 0);
    fraction scale;
    fraction difference;

    for (size_t i = 0; i < dimension_; ++i) {
        for (size_t j = 0; j < dimension_; ++j) {
            fraction::sub(difference, game_frc_(i, j), translation);
            if (difference.sgn() < 0) {
                fmpq_neg(difference.data(), difference.data());
            }
            if (difference > scale) {
                scale = difference;
            }
        }
    }
    if (scale.is_zero()) {
        return false;
    }

    game_dbl_ = linalg::matrix_dbl(dimension_, dimension_);
    fraction normalized;
    for (size_t i = 0; i < dimension_; ++i) {
        for (size_t j = 0; j < dimension_; ++j) {
            fraction::sub(difference, game_frc_(i, j), translation);
            fraction::div(normalized, difference, scale);
            const double converted = normalized.to_dbl();
            if (!std::isfinite(converted) ||
                (converted == 0.0 && !normalized.is_zero()) ||
                std::fpclassify(converted) == FP_SUBNORMAL) {
                return false;
            }
            game_dbl_(i, j) = converted;
        }
    }
    return true;
}

bool find_candidate_unsafe::find(const bitset64& support, size_t support_size)
{
    /*
     * False means the floating result looks clearly incompatible with a
     * candidate. True means "ask the exact solver" for both plausible and
     * numerically suspicious cases. This heuristic can still reject a valid
     * candidate, so exact mode bypasses it completely.
     */
    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    const size_t support_count =
        bs64::extract_set_indices(support, dimension_, support_indices);
    const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
    const size_t non_support_count =
        bs64::extract_set_indices(complement, dimension_, non_support_indices);

    // Build [A'_S -1 | 0; 1^T 0 | 1] directly in reusable scratch.
    const size_t system_dimension = support_size + 1;
    if (linear_system_.rows() != system_dimension) {
        linear_system_ = linalg::matrix_dbl(system_dimension, system_dimension + 1);
    }
    auto& system = linear_system_;
    for (size_t row = 0; row < support_size; ++row) {
        const size_t i = support_indices[row];
        for (size_t column = 0; column < support_size; ++column) {
            system(row, column) = game_dbl_(i, support_indices[column]);
        }
        system(row, support_size) = -1.0;
        system(row, system_dimension) = 0.0;
    }
    for (size_t column = 0; column < support_size; ++column) {
        system(support_size, column) = 1.0;
    }
    system(support_size, support_size) = 0.0;
    system(support_size, system_dimension) = 1.0;

    constexpr double kPivotTolerance = 1e-12;
    double minimum_pivot = std::numeric_limits<double>::infinity();
    for (size_t column = 0; column + 1 < system_dimension; ++column) {
        size_t max_row = column;
        double max_value = std::abs(system(column, column));
        for (size_t row = column + 1; row < system_dimension; ++row) {
            const double value = std::abs(system(row, column));
            if (value > max_value) {
                max_value = value;
                max_row = row;
            }
        }

        // An unreliable pivot provides no useful unsafe decision.
        if (!std::isfinite(max_value) || max_value < kPivotTolerance) return true;
        minimum_pivot = std::min(minimum_pivot, max_value);
        if (max_row != column) system.swap_rows(column, max_row);

        const double pivot = system(column, column);
        for (size_t row = column + 1; row < system_dimension; ++row) {
            const double factor = system(row, column) / pivot;
            for (size_t j = column + 1; j <= system_dimension; ++j) {
                system(row, j) -= factor * system(column, j);
            }
            system(row, column) = 0.0;
        }
    }

    const double last_pivot =
        std::abs(system(system_dimension - 1, system_dimension - 1));
    if (!std::isfinite(last_pivot) || last_pivot < kPivotTolerance) return true;
    minimum_pivot = std::min(minimum_pivot, last_pivot);

    double solution[bs64::kMaxBitsetDimension];
    double solution_scale = 1.0;
    for (size_t i = system_dimension; i-- > 0;) {
        double sum = system(i, system_dimension);
        for (size_t j = i + 1; j < system_dimension; ++j) {
            sum -= system(i, j) * solution[j];
        }
        const double value = sum / system(i, i);
        if (!std::isfinite(value)) return true;
        solution[i] = value;
        solution_scale = std::max(solution_scale, std::abs(value));
    }

    constexpr double kProbabilityTolerance = -1e-10;
    constexpr double kGuardMultiplier = 64.0;
    constexpr double kBinary64Epsilon = std::numeric_limits<double>::epsilon();
    const double base_risk =
        kGuardMultiplier * static_cast<double>(system_dimension) * kBinary64Epsilon;
    const auto is_suspicious = [&](double margin, double decision_scale) noexcept {
        const double risk = base_risk * decision_scale;
        const double separation = minimum_pivot * margin;
        return !std::isfinite(risk) || !std::isfinite(separation) ||
               separation <= 0.0 || risk >= separation;
    };

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
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        double rowsum = 0.0;
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            rowsum += game_dbl_(i, support_indices[j_pos]) * solution[j_pos];
        }
        const double gain = rowsum - payoff;
        if (!std::isfinite(gain)) return true;
        if (gain > proposal_threshold) {
            return is_suspicious(gain, outside_scale);
        }
    }

    return true;
}

} // namespace candidate_search
