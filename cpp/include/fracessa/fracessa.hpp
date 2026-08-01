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
#include <fracessa/find_candidate_exact.hpp>
#include <fracessa/find_candidate_unsafe.hpp>
#include <fracessa/find_candidate_verified.hpp>

enum class analysis_mode {
    verified,
    exact,
    unsafe,
};

analysis_mode parse_analysis_mode(std::string_view name);

/*
 * Runs the complete ESS search for one payoff matrix.
 *
 * For each support S, the analyzer performs up to three stages:
 * 1) the selected numerical procedure may reject an invalid support;
 * 2) exact rational arithmetic constructs and verifies a symmetric Nash
 *    candidate on S, including its extended support;
 * 3) exact positive-definiteness and copositivity tests decide ESS stability.
 *
 * Verified mode returns false only after a rigorous one-sided proof. Unsafe mode uses the historical raw-double heuristic,
 * and exact mode bypasses every preliminary search.
 */
class fracessa
{
public:
    fracessa(const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates = false,
             analysis_mode mode = analysis_mode::verified, bool full_support = false,
             bool with_log = false, std::int64_t matrix_id = -1);
    fracessa(const fracessa&) = delete;
    fracessa& operator=(const fracessa&) = delete;
    fracessa(fracessa&&) = delete;
    fracessa& operator=(fracessa&&) = delete;

    size_t ess_count_ = 0;
    // Populated only when with_candidates is true.
    std::vector<candidate> candidates_;

private:
    // fracessa owns the rational game used by stability and the three floating-point procedures.
    // The exact procedure owns one integer-scaled copy.
    linalg::matrix_frc game_matrix_;
    candidate_search::find_candidate_verified find_candidate_verified_;
    candidate_search::find_candidate_unsafe find_candidate_unsafe_;
    candidate_search::find_candidate_exact find_candidate_exact_;
    linalg::matrix_frc bee_matrix_;

    size_t dimension_;

    bool conf_with_candidates_;
    analysis_mode mode_;
    bool conf_full_support_;
    bool conf_with_log_;

    candidate candidate_;

    std::shared_ptr<spdlog::logger> logger_;

    // Test one support through the selected numerical mode and exact stability.
    // True means an exact equilibrium was found and may prune later supports.
    bool analyze_support(bitset64 support, size_t support_size);

    // Count and optionally output the one exact candidate in candidate_.
    // Circular representatives carry their distinct dihedral-orbit size.
    void finalize_candidate(std::optional<size_t> multiplier);

    // Classify the exact candidate already stored in candidate_.
    void check_stability();
};

#endif // FRACESSA_HPP
