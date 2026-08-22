#include <fracessa/fracessa.hpp>
#include <fracessa/bitset.hpp>
#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/support_generator_sat.hpp>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
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

template<class SupportMask, class Callback>
size_t for_each_dihedral_image(const SupportMask& support, size_t dimension, Callback&& callback)
{
    if constexpr (std::is_same_v<SupportMask, support::bitset>) {
        const SupportMask reflected = support::reflect(support, dimension);
        bool reflection_is_rotation = false;
        size_t count = 0;
        SupportMask current = support;
        do {
            reflection_is_rotation |= current == reflected;
            callback(current);
            ++count;
            current = support::rot_one_right(current, dimension);
        } while (current != support);

        if (!reflection_is_rotation) {
            current = reflected;
            do {
                callback(current);
                ++count;
                current = support::rot_one_right(current, dimension);
            } while (current != reflected);
        }
        return count;
    } else {
        SupportMask reflected = support;
        reflected.reflect();
        bool reflection_is_rotation = false;
        size_t count = 0;
        SupportMask current = support;
        do {
            reflection_is_rotation |= current == reflected;
            callback(current);
            ++count;
            current.rotate_one_right();
        } while (current != support);

        if (!reflection_is_rotation) {
            current.copy_from(reflected);
            do {
                callback(current);
                ++count;
                current.rotate_one_right();
            } while (current != reflected);
        }
        return count;
    }
}

template<class SupportMask, class Callback>
void for_each_dihedral_image_pair(const SupportMask& required, const SupportMask& invaders, size_t dimension,
                                  Callback&& callback)
{
    if constexpr (std::is_same_v<SupportMask, support::bitset>) {
        const SupportMask reflected_required = support::reflect(required, dimension);
        const SupportMask reflected_invaders = support::reflect(invaders, dimension);
        bool reflection_is_rotation = false;
        SupportMask current_required = required;
        SupportMask current_invaders = invaders;
        do {
            reflection_is_rotation |= current_required == reflected_required;
            callback(current_required, current_invaders);
            current_required = support::rot_one_right(current_required, dimension);
            current_invaders = support::rot_one_right(current_invaders, dimension);
        } while (current_required != required);

        if (!reflection_is_rotation) {
            current_required = reflected_required;
            current_invaders = reflected_invaders;
            do {
                callback(current_required, current_invaders);
                current_required = support::rot_one_right(current_required, dimension);
                current_invaders = support::rot_one_right(current_invaders, dimension);
            } while (current_required != reflected_required);
        }
    } else {
        SupportMask reflected_required = required;
        SupportMask reflected_invaders = invaders;
        reflected_required.reflect();
        reflected_invaders.reflect();
        bool reflection_is_rotation = false;
        SupportMask current_required = required;
        SupportMask current_invaders = invaders;
        do {
            reflection_is_rotation |= current_required == reflected_required;
            callback(current_required, current_invaders);
            current_required.rotate_one_right();
            current_invaders.rotate_one_right();
        } while (current_required != required);

        if (!reflection_is_rotation) {
            current_required.copy_from(reflected_required);
            current_invaders.copy_from(reflected_invaders);
            do {
                callback(current_required, current_invaders);
                current_required.rotate_one_right();
                current_invaders.rotate_one_right();
            } while (current_required != reflected_required);
        }
    }
}

} // namespace

template<class SupportMask>
basic_analyzer<SupportMask>::basic_analyzer(search_method method, coposit::parsers::parsed_matrix game,
                                            bool with_candidates, bool full_support, bool with_log, std::int64_t matrix_id)
    : candidate_structure_(make_structure(game.matrix.rows()))
    , ess_structure_(make_structure(game.matrix.rows()))
    , game_(std::move(game))
    , full_support_probe_(game_.matrix.rows())
    , exact_candidate_solver_(game_.matrix, game_.denominator)
    , dimension_(game_.matrix.rows())
    , invaders_(make_empty_support(dimension_))
    , conf_with_candidates_(with_candidates)
    , conf_full_support_(full_support)
    , candidate_(make_candidate(dimension_))
    , logger_()
{
    if constexpr (std::is_same_v<SupportMask, support::bitset>) {
        if (conf_with_candidates_) candidates_.reserve(250 * dimension_);
    } else {
        if (conf_with_candidates_) candidates_.reserve(checked_product(250, dimension_, "Candidate reserve"));
    }

    SupportMask full_support_mask = [&]() {
        if constexpr (std::is_same_v<SupportMask, support::bitset>) return support::set_all_n_bits(dimension_);
        else {
            SupportMask result(dimension_);
            result.set_all();
            return result;
        }
    }();

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

    if (conf_full_support_) {
        if (check_full_support_exact()) {
            finish_log();
            return;
        }
    }

    support::SatSupportGenerator generator(dimension_);
    SupportMask current_support = make_empty_support(dimension_);

    try {
        if (game_.compact_circular) {
            using symmetry_type = std::conditional_t<std::is_same_v<SupportMask, support::bitset>,
                                                     support::CircularAffineSymmetry,
                                                     support::CircularAffineSymmetryMultiword>;
            symmetry_type symmetry(game_.matrix);

            for (size_t support_size = 1; support_size < dimension_; ++support_size) {
                generator.start_cardinality(support_size);
                while (generator.take_first(current_support)) {
                    check_full_support_after_singletons(support_size);
                    log_support_size(support_size);
                    const support_result result = analyze_support(current_support, support_size);
                    if (result == support_result::curvature_ceiling || result == support_result::candidate)
                        full_support_forbidden = true;

                    if (result == support_result::curvature_ceiling) {
                        symmetry.for_each_distinct_bracelet_image(current_support,
                                [&](const auto& image, size_t, bool, size_t) {
                            for_each_dihedral_image(image, dimension_, [&](const auto& orbit_image) {
                                generator.add_upper_cone(orbit_image);
                            });
                        });
                    } else if (result == support_result::candidate) {
                        std::optional<candidate_type> solved_candidate;
                        if (conf_with_candidates_) solved_candidate = candidate_;
                        const SupportMask solved_support = candidate_.support;

                        symmetry.for_each_distinct_bracelet_image(solved_support,
                                [&](const auto& image, size_t multiplier_class, bool reflected, size_t right_shifts) {
                            const size_t multiplier = for_each_dihedral_image(image, dimension_, [&](const auto& orbit_image) {
                                generator.add_upper_cone(orbit_image);
                            });
                            if (solved_candidate) {
                                const size_t previous_candidate_id = candidate_.candidate_id;
                                candidate_ = *solved_candidate;
                                candidate_.candidate_id = previous_candidate_id;
                                if constexpr (std::is_same_v<SupportMask, support::bitset>) {
                                    candidate_.support = image;
                                    candidate_.extended_support = symmetry.image_mask(
                                        solved_candidate->extended_support, multiplier_class, reflected, right_shifts);
                                } else {
                                    candidate_.support.copy_from(image);
                                    symmetry.image_mask(solved_candidate->extended_support, multiplier_class, reflected,
                                                        right_shifts, candidate_.extended_support);
                                }
                                for (size_t position = 0; position < dimension_; ++position) {
                                    const size_t image_position = symmetry.image_position(
                                        position, multiplier_class, reflected, right_shifts);
                                    candidate_.vector[image_position] = solved_candidate->vector[position];
                                }
                            }
                            finalize_candidate(multiplier);
                        });
                    } else if (result == support_result::invader_interval) {
                        symmetry.for_each_distinct_bracelet_image(current_support,
                                [&](const auto& image, size_t multiplier_class, bool reflected, size_t right_shifts) {
                            SupportMask image_invaders = make_empty_support(dimension_);
                            if constexpr (std::is_same_v<SupportMask, support::bitset>) {
                                image_invaders = symmetry.image_mask(invaders_, multiplier_class, reflected, right_shifts);
                            } else {
                                symmetry.image_mask(invaders_, multiplier_class, reflected, right_shifts, image_invaders);
                            }
                            for_each_dihedral_image_pair(image, image_invaders, dimension_,
                                    [&](const auto& orbit_support, const auto& orbit_invaders) {
                                generator.add_invader_interval(orbit_support, orbit_invaders);
                            });
                        });
                    } else {
                        symmetry.for_each_distinct_bracelet_image(current_support,
                                [&](const auto& image, size_t, bool, size_t) {
                            for_each_dihedral_image(image, dimension_, [&](const auto& orbit_image) {
                                generator.add_exact(orbit_image, support_size);
                            });
                        });
                    }
                }
            }

            if (conf_with_candidates_) {
                std::sort(candidates_.begin(), candidates_.end(), [](const candidate_type& left, const candidate_type& right) {
                    return left.support_size < right.support_size ||
                           (left.support_size == right.support_size && left.support < right.support);
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

template<class SupportMask>
typename basic_analyzer<SupportMask>::support_result
basic_analyzer<SupportMask>::analyze_support(const SupportMask& support, size_t support_size)
{
    const bool candidate_found = exact_candidate_solver_.find(support, support_size, candidate_, conf_with_candidates_, invaders_);
    if (!exact_candidate_solver_.reduced_hessian_is_negative_definite()) return support_result::curvature_ceiling;
    if (!candidate_found) {
        if constexpr (std::is_same_v<SupportMask, support::bitset>)
            return invaders_ == 0 ? support_result::rejected : support_result::invader_interval;
        else
            return invaders_.empty() ? support_result::rejected : support_result::invader_interval;
    }

    candidate_.support_size = support_size;
    candidate_.support = support;

    log_candidate();
    check_stability();
    return support_result::candidate;
}

template<class SupportMask>
void basic_analyzer<SupportMask>::finalize_candidate(std::optional<size_t> multiplier)
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

    if (conf_with_candidates_) candidates_.push_back(candidate_);
}

template basic_analyzer<support::bitset>::basic_analyzer(
    search_method, coposit::parsers::parsed_matrix, bool, bool, bool, std::int64_t);
template basic_analyzer<support::bitset_multiword>::basic_analyzer(
    search_method, coposit::parsers::parsed_matrix, bool, bool, bool, std::int64_t);
template basic_analyzer<support::bitset>::support_result
basic_analyzer<support::bitset>::analyze_support(const support::bitset&, size_t);
template basic_analyzer<support::bitset_multiword>::support_result
basic_analyzer<support::bitset_multiword>::analyze_support(const support::bitset_multiword&, size_t);
template void basic_analyzer<support::bitset>::finalize_candidate(std::optional<size_t>);
template void basic_analyzer<support::bitset_multiword>::finalize_candidate(std::optional<size_t>);

} // namespace fracessa
