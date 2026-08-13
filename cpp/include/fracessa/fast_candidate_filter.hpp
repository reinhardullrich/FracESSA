#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include <fracessa/bitset.hpp>
#include <fracessa/bitset_multiword.hpp>
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

// Production floating-point candidate prefilter. Every surviving support is verified by exact arithmetic.
class fast_candidate_filter {
public:
    explicit fast_candidate_filter(size_t dimension);
    fast_candidate_filter(const fast_candidate_filter&) = delete;
    fast_candidate_filter& operator=(const fast_candidate_filter&) = delete;
    fast_candidate_filter(fast_candidate_filter&&) = delete;
    fast_candidate_filter& operator=(fast_candidate_filter&&) = delete;

    // Normalize the parser's denominator-cleared integer game once and convert it to binary64.
    void convert_game_matrix(const numeric::matrix_int& integer_game);
    safe_fallback safe_fallback_reason() const noexcept { return safe_fallback_; }

    // False means heuristic rejection. True means exact arithmetic must decide.
    bool passes(const support::bitset& support, size_t support_size);
    bool passes(const support::bitset_multiword& support, size_t support_size);

private:
    template<class SupportMask, class Index>
    bool passes_from_indices(const SupportMask& support, const Index* support_indices, size_t support_count,
                           double* solution, double* scale_ratios, int* pivots);

    size_t dimension_;
    numeric::matrix_dbl game_dbl_;
    numeric::matrix_dbl reduced_system_;
    std::array<double, support::kMaxBitsetDimension> game_scales_small_{};
    std::vector<double> game_scales_large_;
    std::vector<size_t> support_indices_large_;
    std::vector<double> solution_large_;
    std::vector<double> scale_ratios_large_;
    std::vector<int> pivots_large_;
    double* game_scales_;
    safe_fallback safe_fallback_ = safe_fallback::none;
};

} // namespace fracessa::search
