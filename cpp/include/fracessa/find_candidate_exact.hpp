#pragma once

#include <flint/fmpz_mat.h>

#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <fracessa/candidate.hpp>
#include <linalg/flint_style_fraction_free_ldlt.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

class find_candidate_exact {
public:
    explicit find_candidate_exact(const linalg::matrix_frc& game_matrix);
    find_candidate_exact(const find_candidate_exact&) = delete;
    find_candidate_exact& operator=(const find_candidate_exact&) = delete;
    ~find_candidate_exact();

    // True means an exact candidate was found and written to result. The dense vector is optional because stability does not use it.
    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector);

    // Valid after find() returns true. The reduced Hessian is Z^T*A_S*Z, where the columns of Z span the zero-sum directions on the support.
    // Its negative definiteness is the support-only second-order ESS condition.
    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

private:
    fmpz* reduced_entry(size_t row, size_t column) noexcept;
    fmpz* right_hand_side_entry(size_t row) noexcept;
    fmpz* solution_entry(size_t row) noexcept;
    const fmpz* solution_entry(size_t row) const noexcept;
    const fmpz* game_entry(size_t row, size_t column) const noexcept;
    void resize_reduced_system(size_t reduced_dimension);
    void build_reduced_system(const uint8_t* support_indices, size_t reduced_dimension);
    void calculate_integer_payoff(fmpz* value, size_t strategy, size_t reference, const uint8_t* support_indices, size_t reduced_dimension);
    void ensure_candidate_vector(candidate& result) const;

    size_t dimension_;
    size_t reduced_dimension_ = 0;
    bool reduced_hessian_is_negative_definite_ = false;
    linalg::fraction_free_ldlt_workspace ffldlt_workspace_;
    fmpz_mat_t integer_game_;
    fmpz_t game_denominator_;
    fmpz_mat_t reduced_system_;
    fmpz_mat_t right_hand_side_;
    fmpz_mat_t solution_numerators_;
    fmpz_t solution_denominator_;
    fmpz_t reference_numerator_;
    fmpz_t payoff_numerator_;
    fmpz_t payoff_denominator_;
    fmpz_t outside_payoff_numerator_;
};

} // namespace candidate_search
