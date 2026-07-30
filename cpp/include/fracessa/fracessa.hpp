#ifndef FRACESSA_HPP
#define FRACESSA_HPP

#include <vector>
#include <memory>
#include <optional>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <linalg/matrix_fraction.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/matrix_server.hpp>

/*
 * Runs the complete ESS search for one payoff matrix.
 *
 * For each support S, the analyzer performs up to three stages:
 * 1) the optional double filter may reject an apparently invalid support;
 * 2) exact rational arithmetic constructs and verifies a symmetric Nash
 *    candidate on S, including its extended support;
 * 3) exact positive-definiteness and copositivity tests decide ESS stability.
 *
 * Stage 1 is an unsafe speed heuristic, not a mathematical certificate. A
 * positive answer from it means only "continue with the exact test." The
 * constructor runs the search synchronously and leaves the results below.
 */
class fracessa
{
public:
    fracessa(const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false, int matrix_id = -1);

    size_t ess_count_ = 0;
    // Populated only when with_candidates is true.
    std::vector<candidate> candidates_;

private:
    MatrixServer matrix_server_;

    size_t dimension_;

    bool conf_with_candidates_;
    bool conf_exact_;
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

    // The double routine is only an unsafe rejection filter. The rational
    // routine is the authoritative equilibrium test and fills candidate_.
    bool find_candidate_dbl(const bitset64& support, size_t support_size);
    bool find_candidate_frc(const bitset64& support, size_t support_size);
    
    // Classify the exact candidate already stored in candidate_.
    void check_stability();
};

#endif // FRACESSA_HPP
