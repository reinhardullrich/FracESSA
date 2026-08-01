#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/supports.hpp>
#include <linalg/matrix_fraction.hpp>

#include <stdexcept>
#include <string>

/*
 * Core search orchestration.
 *
 * The constructor performs the full enumeration process. Supports are scanned
 * by increasing size, and each support runs through candidate detection plus
 * stability classification. Confirmed ESS supports trigger superset pruning.
 */

analysis_mode parse_analysis_mode(std::string_view name)
{
    if (name == "verified") return analysis_mode::verified;
    if (name == "exact") return analysis_mode::exact;
    if (name == "unsafe") return analysis_mode::unsafe;
    throw std::invalid_argument("Unknown analysis mode '" + std::string(name) + "'; expected verified, exact, or unsafe");
}

fracessa::fracessa(const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates,
                   analysis_mode mode, bool full_support, bool with_log,
                   std::int64_t matrix_id)
    : game_matrix_(matrix)
    , find_candidate_verified_(game_matrix_)
    , find_candidate_unsafe_(game_matrix_)
    , find_candidate_exact_(game_matrix_)
    , dimension_(matrix.rows())
    , conf_with_candidates_(with_candidates)
    , mode_(mode)
    , conf_full_support_(full_support)
    , conf_with_log_(with_log)
    , candidate_()
    , logger_()
{
    if (mode_ == analysis_mode::verified) {
        if (const char* reason = candidate_search::unavailable_reason()) {
            throw std::runtime_error(
                std::string("Verified candidate search is unavailable: ") + reason +
                ". Run with --mode exact for correct exact analysis, or --mode unsafe for heuristic rejection that may miss "
                "candidates or ESS results.");
        }
    }

    if (conf_with_candidates_)
        candidates_.reserve(250 * dimension_);

    if (conf_with_log_) {
        auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>("log/fracessa.log", 20 * 1024 * 1024, 5);
        logger_ = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
        logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
        logger_->set_level(spdlog::level::info);

        logger_->info("");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");
        logger_->info("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");
        if (matrix_id >= 0)
            logger_->info("matrix_id={}", matrix_id);

        logger_->info("n={}", dimension_);
        logger_->info("game matrix:\n{}", game_matrix_.to_log_string());
    }

    if (mode_ == analysis_mode::unsafe) {
        find_candidate_unsafe_.convert_game_matrix();
        if (find_candidate_unsafe_.input_warnings()) mode_ = analysis_mode::exact;
    }

    if (conf_full_support_) {
        const bitset64 full_support_mask = bs64::set_all_n_bits(dimension_);
        if (analyze_support(full_support_mask, dimension_))
            finalize_candidate(is_cs ? std::optional<size_t>{1} : std::nullopt);
        if (ess_count_ > 0)
            return;
    }

    if (is_cs) {
        CircularSupportGenerator generator(dimension_);
        size_t last_logged_support_size = 0;
        // This lambda is the callback. The generator calls it once for each circular representative and supplies both the
        // mask and its size.
        generator.generate([&](bitset64 support, size_t support_size) {
            if (conf_full_support_ && support_size == dimension_)
                return;
            if (conf_with_log_ && support_size != last_logged_support_size) {
                logger_->info("Searching support size {}:", support_size);
                last_logged_support_size = support_size;
            }
            if (analyze_support(support, support_size)) {
                const size_t multiplier = generator.add_forbidden(candidate_.support);
                finalize_candidate(multiplier);
            }
        });
    } else {
        NonCircularSupportGenerator generator(dimension_);
        size_t last_logged_support_size = 0;
        // This lambda is the callback. The generator calls it once for each support and supplies both the mask and its size.
        generator.generate([&](bitset64 support, size_t support_size) {
            if (conf_full_support_ && support_size == dimension_)
                return;
            if (conf_with_log_ && support_size != last_logged_support_size) {
                logger_->info("Searching support size {}:", support_size);
                last_logged_support_size = support_size;
            }
            if (analyze_support(support, support_size)) {
                generator.add_forbidden(candidate_.support);
                finalize_candidate(std::nullopt);
            }
        });
    }
}

bool fracessa::analyze_support(bitset64 support, size_t support_size) {
    switch (mode_) {
    case analysis_mode::verified:
        if (!find_candidate_verified_.find(support, support_size)) return false;
        break;
    case analysis_mode::unsafe:
        if (!find_candidate_unsafe_.find(support, support_size)) return false;
        break;
    case analysis_mode::exact:
        break;
    }
    const bool needs_candidate_vector = conf_with_candidates_ || conf_with_log_;
    if (!find_candidate_exact_.find(
            support, support_size, candidate_, needs_candidate_vector))
        return false;

    candidate_.support_size = support_size;
    candidate_.support = support;

    if (conf_with_log_)
        logger_->info("Found candidate! Check stability:");

    check_stability();
    return true;
}

void fracessa::finalize_candidate(std::optional<size_t> multiplier) {
    candidate_.multiplier = multiplier;
    ++candidate_.candidate_id;

    if (candidate_.is_ess)
        ess_count_ += multiplier.value_or(1);

    if (conf_with_candidates_)
        candidates_.push_back(candidate_);
    if (conf_with_log_) {
        logger_->info("{}", candidate::header());
        logger_->info("{}", candidate_.to_string());
    }
}
