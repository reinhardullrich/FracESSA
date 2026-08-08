#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include <fracessa/bitset64.hpp>
#include <fracessa/bitset_multiword.hpp>
#include <fracessa/exact_candidate_solver.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

/*
 * Independent candidate prefilter for numerical experiments. Changes belong here until measurements justify moving them into
 * production fast filtering.
 */
class test_candidate_filter {
public:
    explicit test_candidate_filter(const linalg::matrix_frc& game_matrix);
    test_candidate_filter(const test_candidate_filter&) = delete;
    test_candidate_filter& operator=(const test_candidate_filter&) = delete;
    test_candidate_filter(test_candidate_filter&&) = delete;
    test_candidate_filter& operator=(test_candidate_filter&&) = delete;

    // Reuse the exact solver's denominator-cleared integer game, normalize it once, and convert it to binary64.
    void convert_game_matrix(const exact_candidate_solver& exact_solver);
    safe_fallback safe_fallback_reason() const noexcept { return safe_fallback_; }

    // False means heuristic rejection. True means exact arithmetic must decide.
    bool passes(const bitset64& support, size_t support_size);
    bool passes(const bitset_multiword& support, size_t support_size);

private:
    template<class SupportMask, class Index>
    bool passes_from_indices(const SupportMask& support, const Index* support_indices, size_t support_count,
                           double* solution, double* scale_ratios, int* pivots);

    size_t dimension_;
    linalg::matrix_dbl game_dbl_;
    linalg::matrix_dbl reduced_system_;
    std::array<double, bs64::kMaxBitsetDimension> game_scales_small_{};
    std::vector<double> game_scales_large_;
    std::vector<size_t> support_indices_large_;
    std::vector<double> solution_large_;
    std::vector<double> scale_ratios_large_;
    std::vector<int> pivots_large_;
    double* game_scales_;
    safe_fallback safe_fallback_ = safe_fallback::none;
};

} // namespace candidate_search
