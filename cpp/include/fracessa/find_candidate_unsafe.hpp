#pragma once

#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

class find_candidate_unsafe {
public:
    explicit find_candidate_unsafe(const linalg::matrix_frc& game_matrix) noexcept;

    // Prepare the normalized binary64 game. False selects exact mode instead.
    bool normalize_game_matrix();

    // False means heuristic rejection. True means exact arithmetic must decide.
    bool find(const bitset64& support, size_t support_size);

private:
    const linalg::matrix_frc& game_frc_;
    size_t dimension_;
    linalg::matrix_dbl game_dbl_;
    linalg::matrix_dbl linear_system_;
};

} // namespace candidate_search
