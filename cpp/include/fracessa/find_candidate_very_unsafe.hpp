#pragma once

#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

struct very_unsafe_input_warnings {
    bool difference_below_pivot_cutoff = false;
    bool distinct_values_collapsed = false;
    bool nonzero_became_zero = false;
    bool non_finite_value = false;
    bool subnormal_value = false;
    bool binary_exponent_range_exceeds_precision = false;

    bool any() const noexcept
    {
        return difference_below_pivot_cutoff || distinct_values_collapsed || nonzero_became_zero || non_finite_value ||
               subnormal_value || binary_exponent_range_exceeds_precision;
    }
};

class find_candidate_very_unsafe {
public:
    explicit find_candidate_very_unsafe(const linalg::matrix_frc& game_matrix) noexcept;

    // Convert the exact game directly to binary64 without affine normalization.
    void convert_game_matrix();

    const very_unsafe_input_warnings& input_warnings() const noexcept { return input_warnings_; }

    // False means heuristic rejection. True means exact arithmetic must decide.
    bool find(const bitset64& support, size_t support_size);

private:
    const linalg::matrix_frc& game_frc_;
    size_t dimension_;
    linalg::matrix_dbl game_dbl_;
    linalg::matrix_dbl linear_system_;
    very_unsafe_input_warnings input_warnings_;
};

} // namespace candidate_search
