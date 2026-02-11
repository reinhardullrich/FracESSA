#ifndef FRACESSA_HPP
#define FRACESSA_HPP

#include <vector>
#include <string>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <linalg/matrix_fraction.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/supports.hpp>
#include <fracessa/matrix_server.hpp>

/*
 * fracessa: high-level orchestrator for the ESS search pipeline.
 *
 * Execution stages per support:
 * 1) fast double filter (necessary inequalities),
 * 2) exact rational solve/check,
 * 3) stability classification (positive definiteness / copositivity criteria).
 *
 * The class also owns candidate logging/output state and support-pruning state.
 */
class fracessa
{
public:
    fracessa(const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false, int matrix_id = -1);

    size_t ess_count_ = 0;
    std::vector<candidate> candidates_;

private:
    MatrixServer matrix_server_;

    size_t dimension_;
    bool is_cs_;

    bool conf_with_candidates_;
    bool conf_exact_;
    bool conf_full_support_;
    bool conf_with_log_;

    candidate candidate_;
    Supports supports_;

    std::shared_ptr<spdlog::logger> logger_;

    void search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime = false);
    
    // No more templates
    bool find_candidate_dbl(const bitset64& support, size_t support_size);
    bool find_candidate_frc(const bitset64& support, size_t support_size);
    
    void check_stability();
};

#endif // FRACESSA_HPP
