#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include <fracessa/support.hpp>
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

/** Potentially incomplete binary64 candidate prefilter. */
class fast_candidate_filter {
public:
    explicit fast_candidate_filter(support::SupportContext& support_context);
    fast_candidate_filter(const fast_candidate_filter&) = delete;
    fast_candidate_filter& operator=(const fast_candidate_filter&) = delete;
    fast_candidate_filter(fast_candidate_filter&&) = delete;
    fast_candidate_filter& operator=(fast_candidate_filter&&) = delete;

    void convert_game_matrix(const numeric::matrix_int& integer_game);
    safe_fallback safe_fallback_reason() const noexcept { return safe_fallback_; }

    bool passes(const support::Support& support, size_t support_size);

private:
    bool passes_from_indices(const size_t* support_indices, size_t support_count,
                             const size_t* outside_indices, size_t outside_count);

    support::SupportContext& support_context_;
    numeric::matrix_dbl game_dbl_;
    numeric::matrix_dbl reduced_system_;
    std::vector<double> game_scales_;
    std::vector<double> beta_;
    std::vector<size_t> active_coordinates_;
    std::vector<size_t> support_indices_;
    std::vector<size_t> outside_indices_;
    std::vector<double> solution_;
    std::vector<double> scale_ratios_;
    std::vector<int> pivots_;
    safe_fallback safe_fallback_ = safe_fallback::none;
};

} // namespace fracessa::search
