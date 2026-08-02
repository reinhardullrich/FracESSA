#ifndef FRACESSA_HPP
#define FRACESSA_HPP

#include <cstdint>
#include <vector>
#include <memory>
#include <optional>
#include <string_view>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <linalg/matrix_fraction.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/find_candidate_safe.hpp>
#include <fracessa/find_candidate_fast.hpp>
#include <fracessa/find_candidate_test.hpp>

enum class search_method {
    fast,
    safe,
    test,
};

search_method parse_search_method(std::string_view name);

/*
 * Runs the complete ESS search for one payoff matrix.
 *
 * For each support S, the analyzer performs up to three stages:
 * 1) fast or experimental test search may heuristically reject an invalid support;
 * 2) exact rational arithmetic constructs and verifies a symmetric Nash
 *    candidate on S, including its extended support;
 * 3) exact positive-definiteness and copositivity tests decide ESS stability.
 *
 * Safe search always starts with exact arithmetic. Fast search removes the game's common denominator, checks its exact integer
 * precision span, equilibrates the complete binary64 game, and solves each reduced symmetric candidate system with Bunch-Kaufman
 * LDL^T. A large span or inconclusive pivot falls back to exact arithmetic. Test search is an independent copy for experiments.
 */
class fracessa
{
public:
    fracessa(search_method method, const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates = false,
             bool full_support = false,
             bool with_log = false, std::int64_t matrix_id = -1);
    fracessa(const fracessa&) = delete;
    fracessa& operator=(const fracessa&) = delete;
    fracessa(fracessa&&) = delete;
    fracessa& operator=(fracessa&&) = delete;

    size_t ess_count_ = 0;
    candidate_search::safe_fallback safe_fallback_ = candidate_search::safe_fallback::none;
    // Populated only when with_candidates is true.
    std::vector<candidate> candidates_;

private:
    // fracessa owns the rational game used by stability. Safe owns its integer copy; fast and test own converted double copies.
    linalg::matrix_frc game_matrix_;
    candidate_search::find_candidate_fast find_candidate_fast_;
    candidate_search::find_candidate_safe find_candidate_safe_;
    candidate_search::find_candidate_test find_candidate_test_;
    linalg::matrix_frc bee_matrix_;

    size_t dimension_;

    bool conf_with_candidates_;
    search_method method_;
    bool conf_full_support_;
    bool conf_with_log_;

    candidate candidate_;

    std::shared_ptr<spdlog::logger> logger_;

    // Test one support through the selected search method and exact stability.
    // True means an exact equilibrium was found and may prune later supports.
    bool analyze_support(bitset64 support, size_t support_size);

    // Count and optionally output the one exact candidate in candidate_.
    // Circular representatives carry their distinct dihedral-orbit size.
    void finalize_candidate(std::optional<size_t> multiplier);

    // Classify the exact candidate already stored in candidate_.
    void check_stability();
};

#endif // FRACESSA_HPP
