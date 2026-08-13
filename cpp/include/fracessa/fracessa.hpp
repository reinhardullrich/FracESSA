#ifndef FRACESSA_HPP
#define FRACESSA_HPP

#include <array>
#include <cstdint>
#include <vector>
#include <memory>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>

#include <coposit/parsers/parsed_matrix.hpp>
#include <linalg/matrix_integer.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/exact_candidate_solver.hpp>
#include <fracessa/fast_candidate_filter.hpp>

/** Candidate-search method selected by the CLI or native caller. */
enum class search_method {
    fast, ///< Potentially incomplete binary64 candidate filter followed by exact candidate verification and stability.
    safe, ///< Complete exact candidate search and exact stability.
};

/**
 * Convert `fast` or `safe` to the corresponding search method.
 *
 * @param name Exact lowercase method name.
 * @return Parsed search method.
 * @throws std::invalid_argument if `name` is not one of the two public names.
 */
search_method parse_search_method(std::string_view name);

namespace spdlog {
class logger;
}

/**
 * Runs the configured ESS search for one payoff matrix.
 *
 * For each support S, the analyzer performs up to three stages:
 * 1) the fast filter may heuristically reject an invalid support;
 * 2) exact integer arithmetic constructs and verifies a symmetric Nash candidate on S, including its extended support; only a
 *    successful candidate's probabilities and payoff are materialized as exact rationals for output;
 * 3) exact reduced-Hessian inertia and copositivity tests decide ESS stability.
 *
 * Safe search starts directly with the exact candidate solver. The fast filter removes the game's common denominator, checks its
 * exact integer precision span, equilibrates the complete binary64 game, and solves each reduced symmetric candidate system with
 * Bunch-Kaufman LDL^T. A large span or inconclusive pivot falls back to exact arithmetic.
 *
 * Construction runs the configured analysis synchronously. The public count and structure fields are always populated. Individual
 * candidate rows are retained only when `with_candidates` is true.
 *
 * @tparam SupportMask `bitset64` through dimension 64 or `bitset_multiword` for larger dimensions.
 */
template<class SupportMask>
class basic_fracessa
{
public:
    using candidate_type = basic_candidate<SupportMask>;
    using structure_type = std::conditional_t<std::is_same_v<SupportMask, bitset64>,
                                              std::array<size_t, bs64::kMaxBitsetDimension + 1>, std::vector<size_t>>;

    /**
     * Analyze one already parsed symmetric payoff matrix.
     *
     * @param method Candidate-search route.
     * @param game Exact integer-scaled square payoff matrix, retained positive denominator, and compact-circular marker. The analyzer
     *        takes ownership without rebuilding the exact representation.
     * @param with_candidates Retain representative candidate rows in `candidates_`; counts are unaffected by this setting.
     * @param full_support Check the full support first. Stop if the selected route finds it and exact stability accepts it as an ESS;
     *        otherwise continue without checking it twice.
     * @param with_log Write a diagnostic trace to `log/fracessa.log`.
     * @param matrix_id Signed identifier written to the diagnostic log.
     */
    basic_fracessa(search_method method, coposit::parsers::parsed_matrix game, bool with_candidates = false,
                   bool full_support = false, bool with_log = false, std::int64_t matrix_id = -1);
    basic_fracessa(const basic_fracessa&) = delete;
    basic_fracessa& operator=(const basic_fracessa&) = delete;
    basic_fracessa(basic_fracessa&&) = delete;
    basic_fracessa& operator=(basic_fracessa&&) = delete;

    /// Candidate count found by the selected method, including circular multipliers.
    size_t candidate_count_ = 0;
    /// ESS count found by the selected method, including circular multipliers.
    size_t ess_count_ = 0;
    /// Candidate counts indexed by support size; index zero is unused.
    structure_type candidate_structure_;
    /// ESS counts indexed by support size; index zero is unused.
    structure_type ess_structure_;
    /// Whole-matrix fast-to-safe fallback reason; per-support exact retries do not set it.
    candidate_search::safe_fallback safe_fallback_ = candidate_search::safe_fallback::none;
    /// Representative candidate rows. Populated only when `with_candidates` is true.
    std::vector<candidate_type> candidates_;

private:
    // Own one parser-produced integer game and denominator. Exact work and symmetry detection reference it directly; the fast filter
    // owns only its converted double copy.
    coposit::parsers::parsed_matrix game_;
    candidate_search::fast_candidate_filter fast_candidate_filter_;
    candidate_search::exact_candidate_solver exact_candidate_solver_;
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
