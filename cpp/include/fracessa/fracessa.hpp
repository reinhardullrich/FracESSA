#ifndef FRACESSA_HPP
#define FRACESSA_HPP

#include <array>
#include <cstdint>
#include <vector>
#include <memory>
#include <optional>
#include <string_view>
#include <type_traits>

#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_integer.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/exact_candidate_solver.hpp>
#include <fracessa/fast_candidate_filter.hpp>
#include <fracessa/test_candidate_filter.hpp>

enum class search_method {
    fast,
    safe,
    test,
};

search_method parse_search_method(std::string_view name);

namespace spdlog {
class logger;
}

/*
 * Runs the complete ESS search for one payoff matrix.
 *
 * For each support S, the analyzer performs up to three stages:
 * 1) the fast or experimental test filter may heuristically reject an invalid support;
 * 2) exact integer arithmetic constructs and verifies a symmetric Nash candidate on S, including its extended support; only a
 *    successful candidate's probabilities and payoff are materialized as exact rationals for output;
 * 3) exact reduced-Hessian inertia and copositivity tests decide ESS stability.
 *
 * Safe search starts directly with the exact candidate solver. The fast filter removes the game's common denominator, checks its
 * exact integer precision span, equilibrates the complete binary64 game, and solves each reduced symmetric candidate system with
 * Bunch-Kaufman LDL^T. A large span or inconclusive pivot falls back to exact arithmetic. Test search is an independent copy for
 * experiments.
 */
template<class SupportMask>
class basic_fracessa
{
public:
    using candidate_type = basic_candidate<SupportMask>;
    using structure_type = std::conditional_t<std::is_same_v<SupportMask, bitset64>,
                                              std::array<size_t, bs64::kMaxBitsetDimension + 1>, std::vector<size_t>>;

    basic_fracessa(search_method method, const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates = false,
                   bool full_support = false, bool with_log = false, std::int64_t matrix_id = -1);
    basic_fracessa(const basic_fracessa&) = delete;
    basic_fracessa& operator=(const basic_fracessa&) = delete;
    basic_fracessa(basic_fracessa&&) = delete;
    basic_fracessa& operator=(basic_fracessa&&) = delete;

    size_t candidate_count_ = 0;
    size_t ess_count_ = 0;
    structure_type candidate_structure_;
    structure_type ess_structure_;
    candidate_search::safe_fallback safe_fallback_ = candidate_search::safe_fallback::none;
    // Populated only when with_candidates is true.
    std::vector<candidate_type> candidates_;

private:
    // Keep the parsed rational game for logging and detection of extra circular symmetries. The exact solver owns the denominator-cleared
    // integer copy used by candidate and stability work; the fast and test filters own converted double copies.
    linalg::matrix_frc game_matrix_;
    candidate_search::fast_candidate_filter fast_candidate_filter_;
    candidate_search::exact_candidate_solver exact_candidate_solver_;
    candidate_search::test_candidate_filter test_candidate_filter_;
    linalg::matrix_int scaled_reduced_b_;

    size_t dimension_;

    bool conf_with_candidates_;
    search_method method_;
    bool conf_full_support_;

    candidate_type candidate_;

    std::shared_ptr<spdlog::logger> logger_;
    size_t last_logged_support_size_ = 0;

    // Keep logging calls semantic so the search and stability code do not own formatting details.
    void start_log(search_method requested_method, bool is_cs, std::int64_t matrix_id);
    void log_support_size(size_t support_size);
    void log_candidate();
    void log_reduced_b(const SupportMask& outside_best_replies);
    void set_stability_result(bool is_ess, std::string_view reason);
    void finish_log();

    // Test one support through the selected search method and exact stability.
    // True means an exact equilibrium was found and may prune later supports.
    bool analyze_support(const SupportMask& support, size_t support_size);

    // Count and optionally output the candidate in candidate_. Circular rows carry only their distinct dihedral-orbit size.
    void finalize_candidate(std::optional<size_t> multiplier);

    // Classify the exact candidate already stored in candidate_.
    void check_stability();

    static structure_type make_structure(size_t dimension)
    {
        if constexpr (std::is_same_v<SupportMask, bitset64>) return {};
        else return structure_type(dimension + 1, 0);
    }

    static candidate_type make_candidate(size_t dimension)
    {
        if constexpr (std::is_same_v<SupportMask, bitset64>) return {};
        else return candidate_type(dimension);
    }
};

using fracessa = basic_fracessa<bitset64>;
using multiword_fracessa = basic_fracessa<bitset_multiword>;

#endif // FRACESSA_HPP
