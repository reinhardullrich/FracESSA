#pragma once

#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <fracessa/candidate.hpp>
#include <linalg/flint_style_fraction_free_ldlt.hpp>
#include <linalg/integer.hpp>
#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_integer.hpp>

namespace candidate_search {

class find_candidate_safe {
public:
    explicit find_candidate_safe(const linalg::matrix_frc& game_matrix);
    find_candidate_safe(const find_candidate_safe&) = delete;
    find_candidate_safe& operator=(const find_candidate_safe&) = delete;

    // Test whether the exact integer-scaled game spans at least the requested factor without rebuilding that representation.
    bool precision_span_at_least(unsigned long limit) const;

    // True means an exact candidate was found and written to result. The dense vector is optional because stability does not use it.
    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector);

    // Valid after find() returns true. The reduced Hessian is Z^T*A_S*Z, where the columns of Z span the zero-sum directions on the support.
    // Its negative definiteness is the support-only second-order ESS condition.
    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

private:
    void resize_reduced_system(size_t reduced_dimension);
    void build_reduced_system(const uint8_t* support_indices, size_t reduced_dimension);
    void calculate_integer_payoff(linalg::integer& value, size_t strategy, size_t reference, const uint8_t* support_indices,
                                  size_t reduced_dimension);
    void ensure_candidate_vector(candidate& result) const;

    size_t dimension_;
    size_t reduced_dimension_ = 0;
    bool reduced_hessian_is_negative_definite_ = false;
    linalg::fraction_free_ldlt_workspace ffldlt_workspace_;
    linalg::matrix_int integer_game_;
    linalg::integer game_denominator_;
    linalg::matrix_int reduced_system_;
    linalg::matrix_int right_hand_side_;
    linalg::matrix_int solution_numerators_;
    linalg::integer solution_denominator_;
    linalg::integer reference_numerator_;
    linalg::integer payoff_numerator_;
    linalg::integer payoff_denominator_;
    linalg::integer outside_payoff_numerator_;
};

} // namespace candidate_search
