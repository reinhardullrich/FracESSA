#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/supports.hpp>
#include <linalg/matrix_fraction.hpp>

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>

/*
 * Core search orchestration.
 *
 * The constructor performs the full enumeration process. Supports are scanned
 * by increasing size, and each support runs through candidate detection plus
 * stability classification. Every exact equilibrium support triggers superset
 * pruning; stability determines only whether that candidate is an ESS.
 */

search_method parse_search_method(std::string_view name)
{
    if (name == "fast") return search_method::fast;
    if (name == "safe") return search_method::safe;
    if (name == "test") return search_method::test;
    throw std::invalid_argument("Unknown search method '" + std::string(name) + "'; expected fast, safe, or test");
}

fracessa::fracessa(search_method method, const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates,
                   bool full_support, bool with_log, std::int64_t matrix_id)
    : game_matrix_(matrix)
    , find_candidate_fast_(game_matrix_)
    , find_candidate_safe_(game_matrix_)
    , find_candidate_test_(game_matrix_)
    , dimension_(matrix.rows())
    , conf_with_candidates_(with_candidates)
    , method_(method)
    , conf_full_support_(full_support)
    , conf_with_log_(with_log)
    , candidate_()
    , logger_()
{
    if (conf_with_candidates_ || conf_with_log_)
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

    // Candidate rows are logged only after their final ordering and IDs are known. Logging-only runs retain them temporarily.
    const auto finish_candidate_output = [&]() {
        if (conf_with_log_ && !candidates_.empty()) {
            logger_->info("{}", candidate::header());
            for (const candidate& row : candidates_) logger_->info("{}", row.to_string());
        }
        if (!conf_with_candidates_) candidates_.clear();
    };

    if (method_ == search_method::fast) {
        find_candidate_fast_.convert_game_matrix(find_candidate_safe_);
        safe_fallback_ = find_candidate_fast_.safe_fallback_reason();
        if (safe_fallback_ != candidate_search::safe_fallback::none) method_ = search_method::safe;
    } else if (method_ == search_method::test) {
        find_candidate_test_.convert_game_matrix(find_candidate_safe_);
        safe_fallback_ = find_candidate_test_.safe_fallback_reason();
        if (safe_fallback_ != candidate_search::safe_fallback::none) method_ = search_method::safe;
    }

    if (conf_full_support_) {
        const bitset64 full_support_mask = bs64::set_all_n_bits(dimension_);
        if (analyze_support(full_support_mask, dimension_))
            finalize_candidate(is_cs ? std::optional<size_t>{1} : std::nullopt);
        if (ess_count_ > 0) {
            finish_candidate_output();
            return;
        }
    }

    if (is_cs) {
        CircularSupportGeneratorV3 generator(dimension_);
        std::optional<CircularAffineSymmetry> symmetry(std::in_place, game_matrix_);
        if (symmetry->multiplier_class_count() == 1) symmetry.reset();
        size_t last_logged_support_size = 0;
        // This lambda is the callback. The generator calls it once for each circular representative and supplies both the
        // mask and its size.
        generator.generate([&](bitset64 support, size_t support_size) {
            if (conf_full_support_ && support_size == dimension_)
                return;
            if (symmetry && !symmetry->is_representative(support))
                return;
            if (conf_with_log_ && support_size != last_logged_support_size) {
                logger_->info("Searching support size {}:", support_size);
                last_logged_support_size = support_size;
            }
            if (analyze_support(support, support_size)) {
                if (symmetry) {
                    std::optional<candidate> solved_candidate;
                    if (conf_with_candidates_ || conf_with_log_) solved_candidate = candidate_;

                    symmetry->for_each_distinct_bracelet_image(candidate_.support,
                            [&](bitset64 image, size_t multiplier_class, bool reflected, size_t right_shifts) {
                        const size_t multiplier = generator.add_forbidden(image);
                        if (solved_candidate) {
                            const size_t previous_candidate_id = candidate_.candidate_id;
                            candidate_ = *solved_candidate;
                            candidate_.candidate_id = previous_candidate_id;
                            candidate_.support = image;
                            candidate_.extended_support = symmetry->image_mask(
                                solved_candidate->extended_support, multiplier_class, reflected, right_shifts);
                            for (size_t position = 0; position < dimension_; ++position) {
                                const size_t image_position = symmetry->image_position(
                                    position, multiplier_class, reflected, right_shifts);
                                candidate_.vector(image_position, 0) = solved_candidate->vector(position, 0);
                            }
                        }
                        finalize_candidate(multiplier);
                    });
                } else {
                    finalize_candidate(generator.add_forbidden(candidate_.support));
                }
            }
        });

        if ((conf_with_candidates_ || conf_with_log_) && symmetry) {
            std::sort(candidates_.begin(), candidates_.end(), [](const candidate& left, const candidate& right) {
                return left.support_size < right.support_size || (left.support_size == right.support_size && left.support < right.support);
            });
            for (size_t index = 0; index < candidates_.size(); ++index) candidates_[index].candidate_id = index + 1;
        }
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

    finish_candidate_output();
}

bool fracessa::analyze_support(bitset64 support, size_t support_size) {
    if (method_ == search_method::fast && !find_candidate_fast_.find(support, support_size)) return false;
    if (method_ == search_method::test && !find_candidate_test_.find(support, support_size)) return false;
    const bool needs_candidate_vector = conf_with_candidates_ || conf_with_log_;
    if (!find_candidate_safe_.find(
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

    if (conf_with_candidates_ || conf_with_log_)
        candidates_.push_back(candidate_);
}
