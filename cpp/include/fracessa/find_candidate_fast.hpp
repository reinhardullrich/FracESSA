#pragma once

#include <array>
#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

class find_candidate_safe;

// Production floating-point candidate prefilter. Every surviving support is verified by exact arithmetic.
class find_candidate_fast {
public:
    explicit find_candidate_fast(const linalg::matrix_frc& game_matrix) noexcept;

    // Reuse the safe solver's denominator-cleared integer game, normalize it once, and convert it to binary64.
    void convert_game_matrix(const find_candidate_safe& safe_search);
    bool requires_safe_fallback() const noexcept { return requires_safe_fallback_; }

    // False means heuristic rejection. True means exact arithmetic must decide.
    bool find(const bitset64& support, size_t support_size);

private:
    size_t dimension_;
    linalg::matrix_dbl game_dbl_;
    linalg::matrix_dbl reduced_system_;
    std::array<double, bs64::kMaxBitsetDimension> game_scales_{};
    bool requires_safe_fallback_ = false;
};

} // namespace candidate_search
