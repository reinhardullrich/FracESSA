#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include <fracessa/bitset.hpp>
#include <fracessa/types.hpp>

namespace fracessa::search {

enum class safe_fallback : std::uint8_t {
    none,
    precision_span,
    equilibration_invalid,
    equilibration_non_convergence,
};

constexpr std::string_view safe_fallback_name(safe_fallback fallback) noexcept
{
    switch (fallback) {
    case safe_fallback::none: return "";
    case safe_fallback::precision_span: return "precision_span";
    case safe_fallback::equilibration_invalid: return "equilibration_invalid";
    case safe_fallback::equilibration_non_convergence: return "equilibration_non_convergence";
    }
    return "";
}

/** Binary64 scheduling probe used once for the full simplex in model A2. */
class fast_candidate_filter {
public:
    explicit fast_candidate_filter(size_t dimension);
    fast_candidate_filter(const fast_candidate_filter&) = delete;
    fast_candidate_filter& operator=(const fast_candidate_filter&) = delete;
    fast_candidate_filter(fast_candidate_filter&&) = delete;
    fast_candidate_filter& operator=(fast_candidate_filter&&) = delete;

    void convert_game_matrix(const numeric::matrix_int& integer_game);
    safe_fallback safe_fallback_reason() const noexcept { return safe_fallback_; }

    /**
     * Preserve A1's exact check when the full support may be an ESS or binary64 is uncertain. Also return true when A2's reliable estimated
     * inertia ceiling is at most ceil(n/2).
     */
    bool full_support_needs_exact_check();

private:
    bool full_support_needs_exact_check(double* solution, double* scale_ratios, int* pivots);

    size_t dimension_;
    numeric::matrix_dbl game_dbl_;
    numeric::matrix_dbl reduced_system_;
    std::array<double, support::kMaxBitsetDimension> game_scales_small_{};
    std::vector<double> game_scales_large_;
    std::vector<double> solution_large_;
    std::vector<double> scale_ratios_large_;
    std::vector<int> pivots_large_;
    double* game_scales_;
    safe_fallback safe_fallback_ = safe_fallback::none;
};

} // namespace fracessa::search
