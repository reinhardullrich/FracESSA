#pragma once

#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

class find_candidate_safe;

/*
 * Independent copy of the fast candidate prefilter for numerical experiments. Changes belong here until measurements justify
 * moving them into production fast search.
 */
class find_candidate_test {
public:
    explicit find_candidate_test(const linalg::matrix_frc& game_matrix) noexcept;

    // Reuse the safe solver's exact integer game to decide whether this matrix must use safe search, then convert to binary64.
    void convert_game_matrix(const find_candidate_safe& safe_search);
    bool requires_safe_fallback() const noexcept { return requires_safe_fallback_; }

    // False means heuristic rejection. True means exact arithmetic must decide.
    bool find(const bitset64& support, size_t support_size);

private:
    const linalg::matrix_frc& game_frc_;
    size_t dimension_;
    linalg::matrix_dbl game_dbl_;
    linalg::matrix_dbl linear_system_;
    bool requires_safe_fallback_ = false;
};

} // namespace candidate_search
