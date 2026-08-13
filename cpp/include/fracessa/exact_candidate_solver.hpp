#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <fracessa/bitset64.hpp>
#include <fracessa/candidate.hpp>
#include <linalg/fraction_free_ldlt_kkt.hpp>
#include <linalg/integer.hpp>
#include <linalg/matrix_integer.hpp>

namespace candidate_search {

class exact_candidate_solver {
public:
    exact_candidate_solver(const linalg::matrix_int& integer_game, const linalg::integer& game_denominator);
    exact_candidate_solver(const exact_candidate_solver&) = delete;
    exact_candidate_solver& operator=(const exact_candidate_solver&) = delete;

    // True means an exact candidate was found and written to result. Public rational output is optional because stability does not use it.
    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_output);
    bool find(const bitset_multiword& support, size_t support_size, multiword_candidate& result, bool materialize_output);

    // Valid after find() returns true. The reduced Hessian is Z^T*A_S*Z, where the columns of Z span the zero-sum directions on the support.
    // Its negative definiteness is the support-only second-order ESS condition.
    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

    // Build an integer matrix that is a positive multiple of Bomze's final reduced B^(r) matrix and return the number of outside
    // best replies.
    // Method: eliminate the unrestricted support block through one exact Schur complement; no inverse is formed.
    // Valid immediately after find() succeeds for the same support and reports a negative-definite reduced Hessian.
    size_t build_scaled_reduced_b(bitset64 support, bitset64 outside_best_replies, linalg::matrix_int& result);
    size_t build_scaled_reduced_b(const bitset_multiword& support, const bitset_multiword& outside_best_replies,
                                  linalg::matrix_int& result);

private:
    void resize_reduced_system(size_t reduced_dimension);
    template<class Index> void build_reduced_system(const Index* support_indices, size_t reduced_dimension);
    linalg::integer::const_reference reduced_entry(size_t reference, size_t row, size_t column);
    template<class Index>
    void calculate_integer_payoff(linalg::integer& value, size_t strategy, size_t reference, const Index* support_indices,
                                  size_t reduced_dimension);
    template<class SupportMask> void ensure_candidate_vector(basic_candidate<SupportMask>& result) const;
    template<class SupportMask, class Index>
    bool find_from_indices(const SupportMask& support, size_t support_size, basic_candidate<SupportMask>& result,
                           bool materialize_output, const Index* support_indices, size_t support_count,
                           const Index* non_support_indices, size_t non_support_count);
    template<class Index>
    size_t build_scaled_reduced_b_from_indices(const Index* support_indices, size_t support_size, const Index* outside_indices,
                                               size_t outside_size, linalg::matrix_int& result);

    size_t dimension_;
    // Dense (reference, row, column) storage makes a cache hit one byte check and one direct array access. Entries are still
    // calculated lazily; the deliberately unused combinations cost less than two MiB even at dimension 64.
    std::vector<linalg::integer> reduced_entry_cache_;
    std::vector<std::uint8_t> reduced_entry_cache_ready_;
    size_t reduced_system_dimension_ = 0;
    bool reduced_hessian_is_negative_definite_ = false;
    linalg::kkt_fraction_free_ldlt_workspace ffldlt_workspace_;
    const linalg::matrix_int& integer_game_;
    const linalg::integer& game_denominator_;
    linalg::matrix_int reduced_system_;
    linalg::matrix_int right_hand_side_;
    linalg::matrix_int solution_numerators_;
    linalg::matrix_int stability_solution_numerators_;
    linalg::integer solution_denominator_;
    linalg::integer reference_numerator_;
    linalg::integer payoff_numerator_;
    linalg::integer payoff_denominator_;
    linalg::integer outside_payoff_numerator_;
    // Allocated only for dimensions above 64. The one-word path keeps its existing stack arrays.
    std::vector<size_t> support_indices_large_;
    std::vector<size_t> non_support_indices_large_;
};

} // namespace candidate_search
