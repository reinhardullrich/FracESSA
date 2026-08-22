#include <fracessa/fracessa.hpp>
#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/support_generator_circular.hpp>
#include <fracessa/support_generator_non_circular.hpp>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

namespace fracessa {

search_method parse_search_method(std::string_view name)
{
    if (name == "fast") return search_method::fast;
    if (name == "safe") return search_method::safe;
    throw std::invalid_argument("Unknown search method '" + std::string(name) + "'; expected fast or safe");
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

analyzer::analyzer(search_method method, coposit::parsers::parsed_matrix game,
                   bool with_candidates, bool full_support, bool with_log, std::int64_t matrix_id)
    : candidate_structure_(game.matrix.rows() + 1, 0)
    , ess_structure_(game.matrix.rows() + 1, 0)
    , game_(std::move(game))
    , dimension_(game_.matrix.rows())
    , support_context_(dimension_)
    , fast_candidate_filter_(support_context_)
    , exact_candidate_solver_(game_.matrix, game_.denominator, support_context_)
    , conf_with_candidates_(with_candidates)
    , method_(method)
    , conf_full_support_(full_support)
    , candidate_(support_context_)
    , outside_best_replies_(support_context_.make())
    , logger_()
{
    if (conf_with_candidates_)
        candidates_.reserve(checked_product(250, dimension_, "Candidate reserve"));

    if (method_ == search_method::fast) {
        fast_candidate_filter_.convert_game_matrix(game_.matrix);
        safe_fallback_ = fast_candidate_filter_.safe_fallback_reason();
        if (safe_fallback_ != search::safe_fallback::none) method_ = search_method::safe;
    }

    if (with_log) start_log(method, game_.compact_circular, matrix_id);

    if (conf_full_support_) {
        support::Support full_support_mask = support_context_.make();
        support_context_.set_all(full_support_mask);
        log_support_size(dimension_);
        if (analyze_support(full_support_mask, dimension_))
            finalize_candidate(game_.compact_circular ? std::optional<size_t>{1} : std::nullopt);
        if (ess_count_ > 0) {
            finish_log();
            return;
        }
    }

    if (game_.compact_circular) {
        support::CircularSupportGenerator generator(support_context_);
        std::optional<support::CircularAffineSymmetry> symmetry;
        symmetry.emplace(game_.matrix, support_context_);
        if (symmetry->multiplier_class_count() == 1) symmetry.reset();

        generator.generate([&](const support::Support& generated_support, size_t support_size) {
            if (conf_full_support_ && support_size == dimension_) return;
            if (symmetry && !symmetry->is_representative(generated_support)) return;
            log_support_size(support_size);
            if (!analyze_support(generated_support, support_size)) return;

            if (symmetry) {
                std::optional<candidate_type> solved_candidate;
                if (conf_with_candidates_) {
                    solved_candidate.emplace(support_context_);
                    solved_candidate->copy_from(candidate_);
                }

                symmetry->for_each_distinct_bracelet_image(
                    candidate_.support,
                    [&](const support::Support& image, size_t multiplier_class, bool reflected, size_t right_shifts) {
                        const size_t multiplier = generator.add_forbidden_orbit(image);
                        if (solved_candidate) {
                            const size_t previous_candidate_id = candidate_.candidate_id;
                            candidate_.copy_from(*solved_candidate);
                            candidate_.candidate_id = previous_candidate_id;
                            support_context_.copy(candidate_.support, image);
                            symmetry->image_mask(solved_candidate->extended_support, multiplier_class, reflected,
                                                 right_shifts, candidate_.extended_support);
                            for (size_t position = 0; position < dimension_; ++position) {
                                const size_t image_position = symmetry->image_position(
                                    position, multiplier_class, reflected, right_shifts);
                                candidate_.vector[image_position] = solved_candidate->vector[position];
                            }
                        }
                        finalize_candidate(multiplier);
                    });
            } else {
                finalize_candidate(generator.add_forbidden_orbit(candidate_.support));
            }
        });

        if (conf_with_candidates_ && symmetry) {
            std::sort(candidates_.begin(), candidates_.end(), [&](const candidate_type& left, const candidate_type& right) {
                return left.support_size < right.support_size
                       || (left.support_size == right.support_size && support_context_.less(left.support, right.support));
            });
            for (size_t index = 0; index < candidates_.size(); ++index) candidates_[index].candidate_id = index + 1;
        }
    } else {
        support::NonCircularSupportGenerator generator(support_context_);
        generator.generate([&](const support::Support& generated_support, size_t support_size) {
            if (conf_full_support_ && support_size == dimension_) return;
            log_support_size(support_size);
            if (analyze_support(generated_support, support_size)) {
                generator.add_forbidden(candidate_.support);
                finalize_candidate(std::nullopt);
            }
        });
    }

    finish_log();
}

bool analyzer::analyze_support(const support::Support& support, size_t support_size)
{
    if (method_ == search_method::fast && !fast_candidate_filter_.passes(support, support_size)) return false;
    if (!exact_candidate_solver_.find(support, support_size, candidate_, conf_with_candidates_)) return false;

    candidate_.support_size = support_size;
    support_context_.copy(candidate_.support, support);

    log_candidate();
    check_stability();
    return true;
}

void analyzer::finalize_candidate(std::optional<size_t> multiplier)
{
    candidate_.multiplier = multiplier;
    add_checked(candidate_.candidate_id, 1, "Candidate ID");

    const size_t count = multiplier.value_or(1);
    add_checked(candidate_count_, count, "Candidate count");
    add_checked(candidate_structure_[candidate_.support_size], count, "Candidate structure");
    if (candidate_.is_ess) {
        add_checked(ess_count_, count, "ESS count");
        add_checked(ess_structure_[candidate_.support_size], count, "ESS structure");
    }

    if (conf_with_candidates_) candidates_.push_back(candidate_.clone());
}

} // namespace fracessa
