#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/circular_support_generator.hpp>
#include <fracessa/non_circular_support_generator.hpp>
#include <linalg/matrix_fraction.hpp>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

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

namespace {

size_t checked_product(size_t left, size_t right, const char* field)
{
    if (left != 0 && right > std::numeric_limits<size_t>::max() / left)
        throw std::overflow_error(std::string(field) + " overflow");
    return left * right;
}

void add_checked(size_t& destination, size_t value, const char* field)
{
    if (value > std::numeric_limits<size_t>::max() - destination)
        throw std::overflow_error(std::string(field) + " overflow");
    destination += value;
}

} // namespace

template<class SupportMask>
basic_fracessa<SupportMask>::basic_fracessa(search_method method, const linalg::matrix_frc& matrix, bool is_cs,
                                            bool with_candidates, bool full_support, bool with_log, std::int64_t matrix_id)
    : candidate_structure_(make_structure(matrix.rows()))
    , ess_structure_(make_structure(matrix.rows()))
    , game_matrix_(matrix)
    , fast_candidate_filter_(game_matrix_)
    , exact_candidate_solver_(game_matrix_)
    , test_candidate_filter_(game_matrix_)
    , dimension_(matrix.rows())
    , conf_with_candidates_(with_candidates)
    , method_(method)
    , conf_full_support_(full_support)
    , conf_with_log_(with_log)
    , candidate_(make_candidate(matrix.rows()))
    , logger_()
{
    if constexpr (std::is_same_v<SupportMask, bitset64>) {
        if (conf_with_candidates_ || conf_with_log_) candidates_.reserve(250 * dimension_);
    } else {
        if (conf_with_candidates_ || conf_with_log_) candidates_.reserve(checked_product(250, dimension_, "Candidate reserve"));
    }

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
            logger_->info("{}", candidate_type::header());
            for (const candidate_type& row : candidates_) logger_->info("{}", row.to_string());
        }
        if (!conf_with_candidates_) candidates_.clear();
    };

    if (method_ == search_method::fast) {
        fast_candidate_filter_.convert_game_matrix(exact_candidate_solver_);
        safe_fallback_ = fast_candidate_filter_.safe_fallback_reason();
        if (safe_fallback_ != candidate_search::safe_fallback::none) method_ = search_method::safe;
    } else if (method_ == search_method::test) {
        test_candidate_filter_.convert_game_matrix(exact_candidate_solver_);
        safe_fallback_ = test_candidate_filter_.safe_fallback_reason();
        if (safe_fallback_ != candidate_search::safe_fallback::none) method_ = search_method::safe;
    }

    if (conf_full_support_) {
        SupportMask full_support_mask = [&]() {
            if constexpr (std::is_same_v<SupportMask, bitset64>) return bs64::set_all_n_bits(dimension_);
            else {
                SupportMask result(dimension_);
                result.set_all();
                return result;
            }
        }();
        if (analyze_support(full_support_mask, dimension_))
            finalize_candidate(is_cs ? std::optional<size_t>{1} : std::nullopt);
        if (ess_count_ > 0) {
            finish_candidate_output();
            return;
        }
    }

    if constexpr (std::is_same_v<SupportMask, bitset64>) {
        if (is_cs) {
            CircularSupportGenerator generator(dimension_);
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
                        std::optional<candidate_type> solved_candidate;
                        if (conf_with_candidates_ || conf_with_log_) solved_candidate = candidate_;

                        symmetry->for_each_distinct_bracelet_image(candidate_.support,
                                [&](bitset64 image, size_t multiplier_class, bool reflected, size_t right_shifts) {
                            const size_t multiplier = generator.add_forbidden_orbit(image);
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
                        finalize_candidate(generator.add_forbidden_orbit(candidate_.support));
                    }
                }
            });

            if ((conf_with_candidates_ || conf_with_log_) && symmetry) {
                std::sort(candidates_.begin(), candidates_.end(), [](const candidate_type& left, const candidate_type& right) {
                    return left.support_size < right.support_size ||
                           (left.support_size == right.support_size && left.support < right.support);
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
    } else {
        if (is_cs) {
            CircularSupportGeneratorMultiword generator(dimension_);
            std::optional<CircularAffineSymmetryMultiword> symmetry(std::in_place, game_matrix_);
            if (symmetry->multiplier_class_count() == 1) symmetry.reset();
            size_t last_logged_support_size = 0;
            generator.generate([&](const bitset_multiword& support, size_t support_size) {
                if (conf_full_support_ && support_size == dimension_) return;
                if (symmetry && !symmetry->is_representative(support)) return;
                if (conf_with_log_ && support_size != last_logged_support_size) {
                    logger_->info("Searching support size {}:", support_size);
                    last_logged_support_size = support_size;
                }
                if (analyze_support(support, support_size)) {
                    if (symmetry) {
                        std::optional<candidate_type> solved_candidate;
                        if (conf_with_candidates_ || conf_with_log_) solved_candidate = candidate_;

                        symmetry->for_each_distinct_bracelet_image(candidate_.support,
                                [&](const bitset_multiword& image, size_t multiplier_class,
                                    bool reflected, size_t right_shifts) {
                            const size_t multiplier = generator.add_forbidden_orbit(image);
                            if (solved_candidate) {
                                const size_t previous_candidate_id = candidate_.candidate_id;
                                candidate_ = *solved_candidate;
                                candidate_.candidate_id = previous_candidate_id;
                                candidate_.support.copy_from(image);
                                symmetry->image_mask(solved_candidate->extended_support, multiplier_class, reflected,
                                                     right_shifts, candidate_.extended_support);
                                for (size_t position = 0; position < dimension_; ++position) {
                                    const size_t image_position = symmetry->image_position(
                                        position, multiplier_class, reflected, right_shifts);
                                    candidate_.vector(image_position, 0) = solved_candidate->vector(position, 0);
                                }
                            }
                            finalize_candidate(multiplier);
                        });
                    } else {
                        finalize_candidate(generator.add_forbidden_orbit(candidate_.support));
                    }
                }
            });

            if ((conf_with_candidates_ || conf_with_log_) && symmetry) {
                std::sort(candidates_.begin(), candidates_.end(), [](const candidate_type& left, const candidate_type& right) {
                    return left.support_size < right.support_size ||
                           (left.support_size == right.support_size && left.support < right.support);
                });
                for (size_t index = 0; index < candidates_.size(); ++index) candidates_[index].candidate_id = index + 1;
            }
        } else {
            NonCircularSupportGeneratorMultiword generator(dimension_);
            size_t last_logged_support_size = 0;
            generator.generate([&](const bitset_multiword& support, size_t support_size) {
                if (conf_full_support_ && support_size == dimension_) return;
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

    finish_candidate_output();
}

template<class SupportMask>
bool basic_fracessa<SupportMask>::analyze_support(const SupportMask& support, size_t support_size) {
    if (method_ == search_method::fast && !fast_candidate_filter_.passes(support, support_size)) return false;
    if (method_ == search_method::test && !test_candidate_filter_.passes(support, support_size)) return false;
    const bool needs_candidate_vector = conf_with_candidates_ || conf_with_log_;
    if (!exact_candidate_solver_.find(
            support, support_size, candidate_, needs_candidate_vector))
        return false;

    candidate_.support_size = support_size;
    candidate_.support = support;

    if (conf_with_log_)
        logger_->info("Found candidate! Check stability:");

    check_stability();
    return true;
}

template<class SupportMask>
void basic_fracessa<SupportMask>::finalize_candidate(std::optional<size_t> multiplier) {
    candidate_.multiplier = multiplier;
    add_checked(candidate_.candidate_id, 1, "Candidate ID");

    const size_t count = multiplier.value_or(1);
    add_checked(candidate_count_, count, "Candidate count");
    add_checked(candidate_structure_[candidate_.support_size], count, "Candidate structure");
    if (candidate_.is_ess) {
        add_checked(ess_count_, count, "ESS count");
        add_checked(ess_structure_[candidate_.support_size], count, "ESS structure");
    }

    if (conf_with_candidates_ || conf_with_log_)
        candidates_.push_back(candidate_);
}

template basic_fracessa<bitset64>::basic_fracessa(search_method, const linalg::matrix_frc&, bool, bool, bool, bool, std::int64_t);
template basic_fracessa<bitset_multiword>::basic_fracessa(
    search_method, const linalg::matrix_frc&, bool, bool, bool, bool, std::int64_t);
template bool basic_fracessa<bitset64>::analyze_support(const bitset64&, size_t);
template bool basic_fracessa<bitset_multiword>::analyze_support(const bitset_multiword&, size_t);
template void basic_fracessa<bitset64>::finalize_candidate(std::optional<size_t>);
template void basic_fracessa<bitset_multiword>::finalize_candidate(std::optional<size_t>);
