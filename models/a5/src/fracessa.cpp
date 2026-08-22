#include <fracessa/fracessa.hpp>
#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/support_generator_sat.hpp>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace fracessa {

search_method parse_search_method(std::string_view name)
{
    if (name == "safe") return search_method::safe;
    throw std::invalid_argument("Unknown search method '" + std::string(name) + "'; expected safe");
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

struct search_complete {};

template<class Callback>
size_t for_each_dihedral_image(support::SupportContext& context, const support::Support& source,
                               support::Support& reflected, support::Support& current, Callback&& callback)
{
    context.copy(reflected, source);
    context.reflect(reflected);
    context.copy(current, source);

    bool reflection_is_rotation = false;
    size_t count = 0;
    do {
        reflection_is_rotation |= context.equal(current, reflected);
        callback(current);
        ++count;
        context.rotate_one_right(current);
    } while (!context.equal(current, source));

    if (!reflection_is_rotation) {
        context.copy(current, reflected);
        do {
            callback(current);
            ++count;
            context.rotate_one_right(current);
        } while (!context.equal(current, reflected));
    }
    return count;
}

template<class Callback>
void for_each_dihedral_image_pair(support::SupportContext& context,
                                  const support::Support& required, const support::Support& invaders,
                                  support::Support& reflected_required, support::Support& reflected_invaders,
                                  support::Support& current_required, support::Support& current_invaders,
                                  Callback&& callback)
{
    context.copy(reflected_required, required);
    context.copy(reflected_invaders, invaders);
    context.reflect(reflected_required);
    context.reflect(reflected_invaders);
    context.copy(current_required, required);
    context.copy(current_invaders, invaders);

    bool reflection_is_rotation = false;
    do {
        reflection_is_rotation |= context.equal(current_required, reflected_required);
        callback(current_required, current_invaders);
        context.rotate_one_right(current_required);
        context.rotate_one_right(current_invaders);
    } while (!context.equal(current_required, required));

    if (!reflection_is_rotation) {
        context.copy(current_required, reflected_required);
        context.copy(current_invaders, reflected_invaders);
        do {
            callback(current_required, current_invaders);
            context.rotate_one_right(current_required);
            context.rotate_one_right(current_invaders);
        } while (!context.equal(current_required, reflected_required));
    }
}

} // namespace

analyzer::analyzer(search_method method, coposit::parsers::parsed_matrix game,
                   bool with_candidates, bool full_support, bool with_log, std::int64_t matrix_id)
    : candidate_structure_(game.matrix.rows() + 1, 0)
    , ess_structure_(game.matrix.rows() + 1, 0)
    , game_(std::move(game))
    , dimension_(game_.matrix.rows())
    , support_context_(dimension_)
    , full_support_probe_(dimension_)
    , exact_candidate_solver_(game_.matrix, game_.denominator, support_context_)
    , conf_with_candidates_(with_candidates)
    , conf_full_support_(full_support)
    , candidate_(support_context_)
    , outside_best_replies_(support_context_.make())
    , invaders_(support_context_.make())
    , logger_()
{
    if (conf_with_candidates_) candidates_.reserve(checked_product(250, dimension_, "Candidate reserve"));

    support::Support full_support_mask = support_context_.make();
    support_context_.set_all(full_support_mask);

    if (with_log) start_log(method, game_.compact_circular, matrix_id);

    bool full_support_checked = false;
    bool full_support_forbidden = false;
    bool full_support_probe_attempted = false;
    size_t maximum_support_size = dimension_;
    const auto check_full_support_exact = [&]() {
        const size_t previous_ess_count = ess_count_;
        full_support_checked = true;
        log_support_size(dimension_);
        const support_result result = analyze_support(full_support_mask, dimension_);
        if (result == support_result::curvature_ceiling)
            maximum_support_size = std::min(dimension_, exact_candidate_solver_.reduced_hessian_negative_inertia() + 1);
        else if (result == support_result::candidate)
            finalize_candidate(game_.compact_circular ? std::optional<size_t>{1} : std::nullopt);
        return ess_count_ > previous_ess_count;
    };
    const auto check_full_support_after_singletons = [&](size_t support_size) {
        if (!full_support_checked && !full_support_probe_attempted && support_size == 2) {
            full_support_probe_attempted = true;
            if (!full_support_forbidden) {
                full_support_probe_.convert_game_matrix(game_.matrix);
                if (full_support_probe_.full_support_needs_exact_check() && check_full_support_exact()) throw search_complete{};
            }
        }
        if (support_size > maximum_support_size) throw search_complete{};
    };

    if (conf_full_support_ && check_full_support_exact()) {
        finish_log();
        return;
    }

    support::SatSupportGenerator generator(support_context_);
    support::Support current_support = support_context_.make();

    try {
        if (game_.compact_circular) {
            support::CircularAffineSymmetry symmetry(game_.matrix, support_context_);
            support::Support reflected = support_context_.make();
            support::Support current = support_context_.make();
            support::Support reflected_required = support_context_.make();
            support::Support reflected_invaders = support_context_.make();
            support::Support current_required = support_context_.make();
            support::Support current_invaders = support_context_.make();
            support::Support image_invaders = support_context_.make();
            support::Support solved_support = support_context_.make();

            for (size_t support_size = 1; support_size < dimension_; ++support_size) {
                generator.start_cardinality(support_size);
                while (generator.take_first(current_support)) {
                    check_full_support_after_singletons(support_size);
                    log_support_size(support_size);
                    const support_result result = analyze_support(current_support, support_size);
                    if (result == support_result::curvature_ceiling || result == support_result::candidate)
                        full_support_forbidden = true;

                    if (result == support_result::curvature_ceiling) {
                        symmetry.for_each_distinct_bracelet_image(
                            current_support, [&](const support::Support& image, size_t, bool, size_t) {
                                for_each_dihedral_image(support_context_, image, reflected, current,
                                    [&](const support::Support& orbit_image) { generator.add_upper_cone(orbit_image); });
                            });
                    } else if (result == support_result::candidate) {
                        std::optional<candidate_type> solved_candidate;
                        if (conf_with_candidates_) {
                            solved_candidate.emplace(support_context_);
                            solved_candidate->copy_from(candidate_);
                        }
                        support_context_.copy(solved_support, candidate_.support);

                        symmetry.for_each_distinct_bracelet_image(
                            solved_support,
                            [&](const support::Support& image, size_t multiplier_class, bool image_reflected, size_t right_shifts) {
                                const size_t multiplier = for_each_dihedral_image(
                                    support_context_, image, reflected, current,
                                    [&](const support::Support& orbit_image) { generator.add_upper_cone(orbit_image); });
                                if (solved_candidate) {
                                    const size_t previous_candidate_id = candidate_.candidate_id;
                                    candidate_.copy_from(*solved_candidate);
                                    candidate_.candidate_id = previous_candidate_id;
                                    support_context_.copy(candidate_.support, image);
                                    symmetry.image_mask(solved_candidate->extended_support, multiplier_class, image_reflected,
                                                        right_shifts, candidate_.extended_support);
                                    for (size_t position = 0; position < dimension_; ++position) {
                                        const size_t image_position = symmetry.image_position(
                                            position, multiplier_class, image_reflected, right_shifts);
                                        candidate_.vector[image_position] = solved_candidate->vector[position];
                                    }
                                }
                                finalize_candidate(multiplier);
                            });
                    } else if (result == support_result::invader_interval) {
                        symmetry.for_each_distinct_bracelet_image(
                            current_support,
                            [&](const support::Support& image, size_t multiplier_class, bool image_reflected, size_t right_shifts) {
                                symmetry.image_mask(invaders_, multiplier_class, image_reflected, right_shifts, image_invaders);
                                for_each_dihedral_image_pair(
                                    support_context_, image, image_invaders, reflected_required, reflected_invaders,
                                    current_required, current_invaders,
                                    [&](const support::Support& orbit_support, const support::Support& orbit_invaders) {
                                        generator.add_invader_interval(orbit_support, orbit_invaders);
                                    });
                            });
                    } else {
                        symmetry.for_each_distinct_bracelet_image(
                            current_support, [&](const support::Support& image, size_t, bool, size_t) {
                                for_each_dihedral_image(support_context_, image, reflected, current,
                                    [&](const support::Support& orbit_image) { generator.add_exact(orbit_image, support_size); });
                            });
                    }
                }
            }

            if (conf_with_candidates_) {
                std::sort(candidates_.begin(), candidates_.end(), [&](const candidate_type& left, const candidate_type& right) {
                    return left.support_size < right.support_size
                           || (left.support_size == right.support_size && support_context_.less(left.support, right.support));
                });
                for (size_t index = 0; index < candidates_.size(); ++index) candidates_[index].candidate_id = index + 1;
            }
        } else {
            for (size_t support_size = 1; support_size < dimension_; ++support_size) {
                generator.start_cardinality(support_size);
                while (generator.take_first(current_support)) {
                    check_full_support_after_singletons(support_size);
                    log_support_size(support_size);
                    const support_result result = analyze_support(current_support, support_size);
                    if (result == support_result::curvature_ceiling || result == support_result::candidate)
                        full_support_forbidden = true;

                    if (result == support_result::curvature_ceiling) {
                        generator.add_upper_cone(current_support);
                    } else if (result == support_result::candidate) {
                        generator.add_upper_cone(candidate_.support);
                        finalize_candidate(std::nullopt);
                    } else if (result == support_result::invader_interval) {
                        generator.add_invader_interval(current_support, invaders_);
                    } else {
                        generator.add_exact(current_support, support_size);
                    }
                }
            }
        }
    } catch (const search_complete&) {}

    // A sideways clause can make one cardinality unsatisfiable while a larger one remains satisfiable, so only the exact support
    // bound or completion of every layer ends SAT enumeration.
    if (!full_support_checked && !full_support_forbidden) check_full_support_exact();

    finish_log();
}

analyzer::support_result analyzer::analyze_support(const support::Support& support, size_t support_size)
{
    const bool candidate_found = exact_candidate_solver_.find(support, support_size, candidate_, conf_with_candidates_, invaders_);
    if (!exact_candidate_solver_.reduced_hessian_is_negative_definite()) return support_result::curvature_ceiling;
    if (!candidate_found)
        return support_context_.empty(invaders_) ? support_result::rejected : support_result::invader_interval;

    candidate_.support_size = support_size;
    support_context_.copy(candidate_.support, support);

    log_candidate();
    check_stability();
    return support_result::candidate;
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
