#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/supports.hpp>
#include <linalg/matrix_fraction.hpp>

fracessa::fracessa(const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates, bool exact,
                   bool full_support, bool with_log, int matrix_id)
    : matrix_server_(matrix)
    , dimension_(matrix.rows())
    , conf_with_candidates_(with_candidates)
    , conf_exact_(exact)
    , conf_full_support_(full_support)
    , conf_with_log_(with_log)
    , candidate_()
    , logger_()
{
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
        logger_->info("game matrix:\n{}", matrix_server_.get_game_matrix_frc().to_log_string());
    }

    if (!conf_exact_ && !matrix_server_.initialize_unsafe_filter())
        conf_exact_ = true;

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
    if (!conf_exact_ && !find_candidate_dbl(support, support_size))
        return false;
    if (!find_candidate_frc(support, support_size))
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
