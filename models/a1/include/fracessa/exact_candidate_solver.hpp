#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <fracessa/candidate.hpp>
#include <fracessa/fraction_free_ldlt_kkt.hpp>
#include <fracessa/support.hpp>
#include <fracessa/types.hpp>

namespace fracessa::search {

/** Exact integer candidate solver for one symmetric payoff matrix. */
class exact_candidate_solver {
public:
    exact_candidate_solver(const numeric::matrix_int& integer_game, const numeric::integer& game_denominator,
                           support::SupportContext& support_context);
    exact_candidate_solver(const exact_candidate_solver&) = delete;
    exact_candidate_solver& operator=(const exact_candidate_solver&) = delete;

    bool find(const support::Support& support, size_t support_size, candidate& result, bool materialize_output);

    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

    size_t build_scaled_reduced_b(const support::Support& support, const support::Support& outside_best_replies,
                                  numeric::matrix_int& result);

private:
    void resize_reduced_system(size_t reduced_dimension);
    void build_reduced_system(const size_t* support_indices, size_t reduced_dimension);
    numeric::integer::const_reference reduced_entry(size_t reference, size_t row, size_t column);
    void calculate_integer_payoff(numeric::integer& value, size_t strategy, size_t reference,
                                  const size_t* support_indices, size_t reduced_dimension);
    void ensure_candidate_vector(candidate& result) const;
    bool find_from_indices(const support::Support& support, size_t support_size, candidate& result, bool materialize_output,
                           const size_t* support_indices, size_t support_count, const size_t* non_support_indices,
                           size_t non_support_count);
    size_t build_scaled_reduced_b_from_indices(const size_t* support_indices, size_t support_size,
                                               const size_t* outside_indices, size_t outside_size,
                                               numeric::matrix_int& result);

    support::SupportContext& support_context_;
    std::vector<numeric::integer> reduced_entry_cache_;
    std::vector<std::uint8_t> reduced_entry_cache_ready_;
    size_t reduced_system_dimension_ = 0;
    bool reduced_hessian_is_negative_definite_ = false;
    numeric::kkt_fraction_free_ldlt_workspace ffldlt_workspace_;
    const numeric::matrix_int& integer_game_;
    const numeric::integer& game_denominator_;
    numeric::matrix_int reduced_system_;
    numeric::matrix_int right_hand_side_;
    numeric::matrix_int solution_numerators_;
    numeric::matrix_int stability_solution_numerators_;
    numeric::integer solution_denominator_;
    numeric::integer reference_numerator_;
    numeric::integer payoff_numerator_;
    numeric::integer payoff_denominator_;
    numeric::integer outside_payoff_numerator_;
    std::vector<size_t> support_indices_;
    std::vector<size_t> non_support_indices_;
};

} // namespace fracessa::search
