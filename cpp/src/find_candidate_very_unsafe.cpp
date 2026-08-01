#include <fracessa/find_candidate_very_unsafe.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace candidate_search {
namespace {

constexpr double kPivotCutoff = 1e-12;

struct converted_entry {
    double value;
    const linalg::fraction* exact;
};

} // namespace

find_candidate_very_unsafe::find_candidate_very_unsafe(
    const linalg::matrix_frc& game_matrix) noexcept
    : game_frc_(game_matrix)
    , dimension_(game_matrix.rows())
{
}

void find_candidate_very_unsafe::convert_game_matrix()
{
    input_warnings_ = {};
    game_dbl_ = linalg::matrix_dbl(dimension_, dimension_);

    // Inspect only the upper triangle because every supported game matrix is symmetric.
    // The checks run once per game, never in the exponential support loop.
    std::vector<converted_entry> finite_entries;
    finite_entries.reserve(dimension_ * (dimension_ + 1) / 2);
    int minimum_exponent = std::numeric_limits<int>::max();
    int maximum_exponent = std::numeric_limits<int>::min();

    for (size_t i = 0; i < dimension_; ++i) {
        for (size_t j = 0; j < dimension_; ++j) {
            const linalg::fraction& exact = game_frc_(i, j);
            const double converted = exact.to_dbl();
            game_dbl_(i, j) = converted;
            if (j < i) continue;

            if (converted == 0.0 && !exact.is_zero()) input_warnings_.nonzero_became_zero = true;
            if (!std::isfinite(converted)) {
                input_warnings_.non_finite_value = true;
                continue;
            }
            if (std::fpclassify(converted) == FP_SUBNORMAL) input_warnings_.subnormal_value = true;
            if (converted != 0.0) {
                const int exponent = std::ilogb(std::abs(converted));
                minimum_exponent = std::min(minimum_exponent, exponent);
                maximum_exponent = std::max(maximum_exponent, exponent);
            }
            finite_entries.push_back({converted, &exact});
        }
    }

    if (minimum_exponent != std::numeric_limits<int>::max() &&
        maximum_exponent - minimum_exponent >= std::numeric_limits<double>::digits) {
        input_warnings_.binary_exponent_range_exceeds_precision = true;
    }

    std::sort(finite_entries.begin(), finite_entries.end(),
              [](const converted_entry& left, const converted_entry& right) { return left.value < right.value; });

    for (size_t first = 0; first < finite_entries.size();) {
        size_t next = first + 1;
        while (next < finite_entries.size() && finite_entries[next].value == finite_entries[first].value) {
            if (!(*finite_entries[next].exact == *finite_entries[first].exact)) input_warnings_.distinct_values_collapsed = true;
            ++next;
        }

        if (next < finite_entries.size()) {
            const double difference = finite_entries[next].value - finite_entries[first].value;
            if (difference > 0.0 && difference < kPivotCutoff) input_warnings_.difference_below_pivot_cutoff = true;
        }
        first = next;
    }
}

bool find_candidate_very_unsafe::find(const bitset64& support, size_t support_size)
{
    /*
     * This is the historical raw-double rejection path from revision 32f61679.
     * It intentionally has no normalization, danger veto, or exact fallback for
     * a small pivot. It can therefore reject valid candidates and ESS results.
     */
    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    const size_t support_count =
        bs64::extract_set_indices(support, dimension_, support_indices);
    const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
    const size_t non_support_count =
        bs64::extract_set_indices(complement, dimension_, non_support_indices);

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

        if (max_value < kPivotCutoff) return false;
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

    if (std::abs(system(system_dimension - 1, system_dimension - 1)) <
        kPivotCutoff) {
        return false;
    }

    double solution[bs64::kMaxBitsetDimension];
    constexpr double kProbabilityTolerance = -1e-10;
    for (size_t i = system_dimension; i-- > 0;) {
        double sum = system(i, system_dimension);
        for (size_t j = i + 1; j < system_dimension; ++j) {
            sum -= system(i, j) * solution[j];
        }
        solution[i] = sum / system(i, i);
        if (i < support_size && solution[i] < kProbabilityTolerance) return false;
    }

    const double threshold =
        solution[support_size] + 1e-4 * static_cast<double>(dimension_);
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        double rowsum = 0.0;
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            rowsum += game_dbl_(i, support_indices[j_pos]) * solution[j_pos];
        }
        if (rowsum > threshold) return false;
    }

    return true;
}

} // namespace candidate_search
