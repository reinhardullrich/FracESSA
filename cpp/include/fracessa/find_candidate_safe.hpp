#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include <fracessa/bitset64.hpp>
#include <fracessa/candidate.hpp>
#include <linalg/flint_style_fraction_free_ldlt.hpp>
#include <linalg/integer.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_integer.hpp>

namespace candidate_search {

enum class safe_fallback : std::uint8_t {
    none,
    precision_span,
    equilibration_invalid,
    equilibration_non_convergence,
};

constexpr std::string_view safe_fallback_name(safe_fallback fallback) noexcept
{
    switch (fallback) {
    case safe_fallback::none: return "";
    case safe_fallback::precision_span: return "precision_span";
    case safe_fallback::equilibration_invalid: return "equilibration_invalid";
    case safe_fallback::equilibration_non_convergence: return "equilibration_non_convergence";
    }
    return "";
}

class find_candidate_safe {
public:
    explicit find_candidate_safe(const linalg::matrix_frc& game_matrix);
    find_candidate_safe(const find_candidate_safe&) = delete;
    find_candidate_safe& operator=(const find_candidate_safe&) = delete;

    // Test whether the exact integer-scaled game spans at least the requested factor without rebuilding that representation.
    bool precision_span_at_least(unsigned long limit) const;

    // Remove the common game denominator, normalize the remaining integer matrix by one power of two, convert one symmetric
    // triangle, and mirror it for double search.
    // Returns none after successful conversion; otherwise identifies why the complete matrix must use safe search.
    safe_fallback prepare_normalized_double_game(unsigned long precision_span_limit, linalg::matrix_dbl& result) const;

    // True means an exact candidate was found and written to result. The dense vector is optional because stability does not use it.
    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector);

    // Valid after find() returns true. The reduced Hessian is Z^T*A_S*Z, where the columns of Z span the zero-sum directions on the support.
    // Its negative definiteness is the support-only second-order ESS condition.
    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

    // Build an integer matrix that is a positive multiple of Bomze's final reduced B^(r) matrix and return its dimension |K|.
    // Method: eliminate the unrestricted support block through one exact Schur complement; no inverse is formed.
    // Valid immediately after find() succeeds for the same support and reports a negative-definite reduced Hessian.
    size_t build_scaled_reduced_b(bitset64 support, bitset64 outside_best_replies, linalg::matrix_frc& result);

private:
    bool precision_span_at_least(unsigned long limit, bool include_game_denominator, linalg::integer& maximum) const;
    void resize_reduced_system(size_t reduced_dimension);
    void build_reduced_system(const uint8_t* support_indices, size_t reduced_dimension);
    linalg::integer::const_reference reduced_entry(size_t reference, size_t row, size_t column);
    void calculate_integer_payoff(linalg::integer& value, size_t strategy, size_t reference, const uint8_t* support_indices,
                                  size_t reduced_dimension);
    void ensure_candidate_vector(candidate& result) const;

    size_t dimension_;
    // Dense (reference, row, column) storage makes a cache hit one byte check and one direct array access. Entries are still
    // calculated lazily; the deliberately unused combinations cost less than two MiB even at dimension 64.
    std::vector<linalg::integer> reduced_entry_cache_;
    std::vector<std::uint8_t> reduced_entry_cache_ready_;
    size_t reduced_system_dimension_ = 0;
    bool reduced_hessian_is_negative_definite_ = false;
    linalg::fraction_free_ldlt_workspace ffldlt_workspace_;
    linalg::matrix_int integer_game_;
    linalg::integer game_denominator_;
    linalg::matrix_int reduced_system_;
    linalg::matrix_int right_hand_side_;
    linalg::matrix_int solution_numerators_;
    linalg::matrix_int stability_solution_numerators_;
    linalg::integer solution_denominator_;
    linalg::integer reference_numerator_;
    linalg::integer payoff_numerator_;
    linalg::integer payoff_denominator_;
    linalg::integer outside_payoff_numerator_;
};

} // namespace candidate_search
