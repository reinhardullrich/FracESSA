#ifndef FRACESSA_HPP
#define FRACESSA_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

#include <coposit/parsers/parsed_matrix.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/exact_candidate_solver.hpp>
#include <fracessa/fast_candidate_filter.hpp>
#include <fracessa/support.hpp>
#include <fracessa/types.hpp>

namespace spdlog {
class logger;
}

namespace fracessa {

/** Candidate-search method selected by the CLI or native caller. */
enum class search_method {
    fast, ///< Compatibility name for A1's exact search with one full-support binary64 probe.
    safe, ///< Compatibility name for the same A1 search and exact decisions.
};

search_method parse_search_method(std::string_view name);

/** Runs the configured ESS search for one payoff matrix. */
class analyzer
{
public:
    using candidate_type = candidate;
    using structure_type = std::vector<size_t>;

    analyzer(search_method method, coposit::parsers::parsed_matrix game, bool with_candidates = false,
             bool full_support = false, bool with_log = false, std::int64_t matrix_id = -1);
    analyzer(const analyzer&) = delete;
    analyzer& operator=(const analyzer&) = delete;
    analyzer(analyzer&&) = delete;
    analyzer& operator=(analyzer&&) = delete;

    const support::SupportContext& support_context() const noexcept { return support_context_; }

    /// Candidate count found by the selected method, including circular multipliers.
    size_t candidate_count_ = 0;
    /// ESS count found by the selected method, including circular multipliers.
    size_t ess_count_ = 0;
    /// Candidate counts indexed by support size; index zero is unused.
    structure_type candidate_structure_;
    /// ESS counts indexed by support size; index zero is unused.
    structure_type ess_structure_;
    /// Whole-matrix fast-to-safe fallback reason; per-support exact retries do not set it.
    search::safe_fallback safe_fallback_ = search::safe_fallback::none;
    /// Representative candidate rows. Populated only when `with_candidates` is true.
    std::vector<candidate_type> candidates_;

private:
    enum class support_result : std::uint8_t {
        rejected,
        candidate,
        curvature_ceiling,
    };

    coposit::parsers::parsed_matrix game_;
    size_t dimension_;
    support::SupportContext support_context_;
    search::fast_candidate_filter fast_candidate_filter_;
    search::exact_candidate_solver exact_candidate_solver_;
    numeric::matrix_int scaled_reduced_b_;

    bool conf_with_candidates_;
    search_method method_;
    bool conf_full_support_;

    candidate_type candidate_;
    support::Support outside_best_replies_;

    std::shared_ptr<spdlog::logger> logger_;
    size_t last_logged_support_size_ = 0;

    void start_log(search_method requested_method, bool is_cs, std::int64_t matrix_id);
    void log_support_size(size_t support_size);
    void log_candidate();
    void log_reduced_b(const support::Support& outside_best_replies);
    void set_stability_result(bool is_ess, std::string_view reason);
    void finish_log();

    support_result analyze_support(const support::Support& support, size_t support_size);
    void finalize_candidate(std::optional<size_t> multiplier);
    void check_stability();
};

} // namespace fracessa

#endif // FRACESSA_HPP
