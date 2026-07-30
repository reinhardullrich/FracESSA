#pragma once

#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <fracessa/candidate.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

class find_candidate_exact {
public:
    explicit find_candidate_exact(const linalg::matrix_frc& game_matrix) noexcept;

    // True means an exact candidate was found and written to result.
    bool find(const bitset64& support, size_t support_size, candidate& result);

private:
    const linalg::matrix_frc& game_frc_;
    size_t dimension_;
    linalg::matrix_frc linear_system_;
};

} // namespace candidate_search
