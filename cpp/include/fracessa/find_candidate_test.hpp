#pragma once

#include <array>
#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <fracessa/find_candidate_safe.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

class find_candidate_safe;

/*
 * Independent candidate prefilter for numerical experiments. Changes belong here until measurements justify moving them into
 * production fast search.
 */
class find_candidate_test {
public:
    explicit find_candidate_test(const linalg::matrix_frc& game_matrix) noexcept;

    // Reuse the safe solver's denominator-cleared integer game, normalize it once, and convert it to binary64.
    void convert_game_matrix(const find_candidate_safe& safe_search);
    safe_fallback safe_fallback_reason() const noexcept { return safe_fallback_; }

    // False means heuristic rejection. True means exact arithmetic must decide.
    bool find(const bitset64& support, size_t support_size);

private:
    size_t dimension_;
    linalg::matrix_dbl game_dbl_;
    linalg::matrix_dbl reduced_system_;
    std::array<double, bs64::kMaxBitsetDimension> game_scales_{};
    safe_fallback safe_fallback_ = safe_fallback::none;
};

} // namespace candidate_search
