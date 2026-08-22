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
#include <fracessa/types.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset.hpp>
#include <fracessa/exact_candidate_solver.hpp>
#include <fracessa/fast_candidate_filter.hpp>

namespace spdlog {
class logger;
}

namespace fracessa {

/** A3 exposes one exact search method. */
enum class search_method {
    safe,
};

/**
 * Accept the exact lowercase method name `safe`.
 *
 * @param name Exact lowercase method name.
 * @return Parsed search method.
 * @throws std::invalid_argument if `name` is not `safe`.
 */
search_method parse_search_method(std::string_view name);

/**
 * Runs the configured ESS search for one payoff matrix.
 *
 * Model A3 checks singleton supports, then uses binary64 only to decide whether an early exact full-support check is worthwhile. All
 * pruning and ESS decisions remain exact. Complete exact full-support inertia bounds every ESS support by `negative_inertia + 1`.
 * A full-game exact row comparison then forbids every support containing a weakly pure-dominated strategy. Exact curvature failures
 * and exact Nash equilibria continue to prune upper cones during ordinary ascending enumeration.
 *
 * Construction runs the configured analysis synchronously. The public count and structure fields are always populated. Individual
 * candidate rows are retained only when `with_candidates` is true.
 *
 * @tparam SupportMask `bitset` through dimension 64 or `bitset_multiword` for larger dimensions.
 */
template<class SupportMask>
class basic_analyzer
{
public:
    using candidate_type = basic_candidate<SupportMask>;
    using structure_type = std::conditional_t<std::is_same_v<SupportMask, support::bitset>,
                                              std::array<size_t, support::kMaxBitsetDimension + 1>, std::vector<size_t>>;

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
    basic_analyzer(search_method method, coposit::parsers::parsed_matrix game, bool with_candidates = false,
                   bool full_support = false, bool with_log = false, std::int64_t matrix_id = -1);
    basic_analyzer(const basic_analyzer&) = delete;
    basic_analyzer& operator=(const basic_analyzer&) = delete;
    basic_analyzer(basic_analyzer&&) = delete;
    basic_analyzer& operator=(basic_analyzer&&) = delete;

    /// Candidate count found by the exact method, including circular multipliers.
    size_t candidate_count_ = 0;
    /// ESS count found by the selected method, including circular multipliers.
    size_t ess_count_ = 0;
    /// Candidate counts indexed by support size; index zero is unused.
    structure_type candidate_structure_;
    /// ESS counts indexed by support size; index zero is unused.
    structure_type ess_structure_;
    /// Always `none`; internal floating-point scheduling failures request exact work and are not public safe fallbacks.
    search::safe_fallback safe_fallback_ = search::safe_fallback::none;
    /// Representative candidate rows. Populated only when `with_candidates` is true.
    std::vector<candidate_type> candidates_;

private:
    enum class support_result : std::uint8_t {
        rejected,
        candidate,
        curvature_ceiling,
    };

    // Own one parser-produced integer game and denominator. Exact work and symmetry detection reference it directly.
    coposit::parsers::parsed_matrix game_;
    search::fast_candidate_filter full_support_probe_;
    search::exact_candidate_solver exact_candidate_solver_;
    numeric::matrix_int scaled_reduced_b_;

    size_t dimension_;

    bool conf_with_candidates_;
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

    // Every generated support is exact. An exact curvature ceiling is distinct from a candidate because both prune supersets, but
    // only the candidate contributes to output.
    support_result analyze_support(const SupportMask& support, size_t support_size);

    // Count and optionally output the candidate in candidate_. Circular rows carry only their distinct dihedral-orbit size.
    void finalize_candidate(std::optional<size_t> multiplier);

    // Classify the exact candidate already stored in candidate_.
    void check_stability();

    static structure_type make_structure(size_t dimension)
    {
        if constexpr (std::is_same_v<SupportMask, support::bitset>) return {};
        else return structure_type(dimension + 1, 0);
    }

    static candidate_type make_candidate(size_t dimension)
    {
        if constexpr (std::is_same_v<SupportMask, support::bitset>) return {};
        else return candidate_type(dimension);
    }
};

using analyzer = basic_analyzer<support::bitset>;
using analyzer_multiword = basic_analyzer<support::bitset_multiword>;

} // namespace fracessa

#endif // FRACESSA_HPP
