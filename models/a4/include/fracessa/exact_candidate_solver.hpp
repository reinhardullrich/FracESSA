#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <fracessa/bitset.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/fraction_free_ldlt_kkt.hpp>
#include <fracessa/types.hpp>

namespace fracessa::search {

/**
 * Exact integer candidate solver for one symmetric payoff matrix.
 *
 * The parser has already cleared the matrix's rational denominators. Each support is therefore solved and verified with
 * fraction-free integer arithmetic. Exact rational probabilities and payoffs are materialized only when candidate rows were
 * requested. The retained reduced-Hessian factorization is reused by the later stability reduction.
 */
class exact_candidate_solver {
public:
    exact_candidate_solver(const numeric::matrix_int& integer_game, const numeric::integer& game_denominator);
    exact_candidate_solver(const exact_candidate_solver&) = delete;
    exact_candidate_solver& operator=(const exact_candidate_solver&) = delete;

    /**
     * Return true after finding an exact candidate.
     *
     * On an outside-payoff rejection, `invaders` contains every outside strategy whose exact payoff is strictly greater than the
     * common support payoff. It is empty after every other result. Rational output is optional because stability uses only exact
     * integers.
     */
    bool find(const support::bitset& support, size_t support_size, candidate& result, bool materialize_output,
              support::bitset& invaders);
    bool find(const support::bitset_multiword& support, size_t support_size, candidate_multiword& result, bool materialize_output,
              support::bitset_multiword& invaders);

    /**
     * Return whether the latest support's reduced Hessian is negative definite.
     *
     * A4 records this result before probability and outside-payoff checks. It is therefore valid after every `find()` call, including
     * a false return. Singular systems and nonsingular systems with positive inertia both return false here.
     */
    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

    /// Exact negative inertia of the latest reduced Hessian, including when that Hessian is singular.
    size_t reduced_hessian_negative_inertia() const noexcept { return reduced_hessian_negative_inertia_; }

    /**
     * Build a positive integer multiple of Bomze's final reduced B^(r) matrix and return its dimension.
     *
     * One exact Schur complement eliminates the unrestricted support block without forming an inverse. Call this immediately after
     * `find()` succeeds for the same support and reports a negative-definite reduced Hessian.
     */
    size_t build_scaled_reduced_b(support::bitset support, support::bitset outside_best_replies, numeric::matrix_int& result);
    size_t build_scaled_reduced_b(const support::bitset_multiword& support, const support::bitset_multiword& outside_best_replies,
                                  numeric::matrix_int& result);

private:
    void resize_reduced_system(size_t reduced_dimension);
    template<class Index> void build_reduced_system(const Index* support_indices, size_t reduced_dimension);
    numeric::integer::const_reference reduced_entry(size_t reference, size_t row, size_t column);
    template<class Index>
    void calculate_integer_payoff(numeric::integer& value, size_t strategy, size_t reference, const Index* support_indices,
                                  size_t reduced_dimension);
    template<class SupportMask> void ensure_candidate_vector(basic_candidate<SupportMask>& result) const;
    template<class SupportMask, class Index>
    bool find_from_indices(const SupportMask& support, size_t support_size, basic_candidate<SupportMask>& result,
                           bool materialize_output, const Index* support_indices, size_t support_count,
                           const Index* non_support_indices, size_t non_support_count, SupportMask& invaders);
    template<class Index>
    size_t build_scaled_reduced_b_from_indices(const Index* support_indices, size_t support_size, const Index* outside_indices,
                                               size_t outside_size, numeric::matrix_int& result);

    size_t dimension_;
    // Dense (reference, row, column) storage makes a cache hit one byte check and one direct array access. Entries are still
    // calculated lazily; the deliberately unused combinations cost less than two MiB even at dimension 64.
    std::vector<numeric::integer> reduced_entry_cache_;
    std::vector<std::uint8_t> reduced_entry_cache_ready_;
    size_t reduced_system_dimension_ = 0;
    bool reduced_hessian_is_negative_definite_ = false;
    size_t reduced_hessian_negative_inertia_ = 0;
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
    // Allocated only for dimensions above 64. The one-word path keeps its existing stack arrays.
    std::vector<size_t> support_indices_large_;
    std::vector<size_t> non_support_indices_large_;
};

} // namespace fracessa::search
