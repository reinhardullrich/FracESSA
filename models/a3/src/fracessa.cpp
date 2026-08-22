#include <fracessa/fracessa.hpp>
#include <fracessa/bitset.hpp>
#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/support_generator_circular.hpp>
#include <fracessa/support_generator_non_circular.hpp>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

/*
 * Core search orchestration.
 *
 * The constructor performs the full enumeration process. Supports are scanned
 * by increasing size, and each support runs through candidate detection plus
 * stability classification. Every exact equilibrium support triggers superset
 * pruning; stability determines only whether that candidate is an ESS.
 */

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

template<class Generator, class Callback>
void generate_until_complete(Generator& generator, Callback&& callback)
{
    try {
        generator.generate(std::forward<Callback>(callback));
    } catch (const search_complete&) {}
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
    std::vector<size_t> weakly_dominated_strategies;

    if (with_log) start_log(method, game_.compact_circular, matrix_id);

    bool full_support_checked = false;
    size_t singleton_supports_seen = 0;
    size_t maximum_support_size = dimension_;
    const auto find_weakly_dominated_strategies = [&]() {
        for (size_t strategy = 0; strategy < dimension_; ++strategy) {
            for (size_t dominating_strategy = 0; dominating_strategy < dimension_; ++dominating_strategy) {
                if (dominating_strategy == strategy) continue;

                size_t opponent = 0;
                while (opponent < dimension_
                       && game_.matrix(dominating_strategy, opponent).compare(game_.matrix(strategy, opponent)) >= 0)
                    ++opponent;
                if (opponent != dimension_) continue;

                weakly_dominated_strategies.push_back(strategy);
                break;
            }
        }
    };
    const auto make_singleton = [&](size_t strategy) {
        if constexpr (std::is_same_v<SupportMask, support::bitset>) return support::set_bit_at_pos(SupportMask{0}, strategy);
        else {
            SupportMask singleton(dimension_);
            singleton.set_bit_at_pos(strategy);
            return singleton;
        }
    };
    const auto finishes_singleton_layer = [&](size_t support_size, size_t expected_singletons) {
        return support_size == 1 && ++singleton_supports_seen == expected_singletons;
    };
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
    const auto check_support_size_limit = [&](size_t support_size) {
        if (support_size > maximum_support_size) throw search_complete{};
    };
    const auto probe_full_support_after_singletons = [&](const auto& generator) {
        if (!generator.has_supports_after_singletons()) throw search_complete{};
        if (!full_support_checked) {
            full_support_probe_.convert_game_matrix(game_.matrix);
            if (full_support_probe_.full_support_needs_exact_check() && check_full_support_exact()) throw search_complete{};
        }
        check_support_size_limit(2);
    };

    if (conf_full_support_) {
        if (check_full_support_exact()) {
            finish_log();
            return;
        }
    }

    if constexpr (std::is_same_v<SupportMask, support::bitset>) {
        if (game_.compact_circular) {
            support::CircularSupportGenerator generator(dimension_);
            std::optional<support::CircularAffineSymmetry> symmetry(std::in_place, game_.matrix);
            if (symmetry->multiplier_class_count() == 1) symmetry.reset();
            // This lambda is the callback. The generator calls it once for each circular representative and supplies both the
            // mask and its size.
            generate_until_complete(generator, [&](support::bitset support, size_t support_size) {
                check_support_size_limit(support_size);
                if (support_size == dimension_) {
                    if (!full_support_checked) check_full_support_exact();
                    return;
                }
                if (symmetry && !symmetry->is_representative(support))
                    return;
                log_support_size(support_size);
                const support_result result = analyze_support(support, support_size);
                if (result == support_result::curvature_ceiling) {
                    if (symmetry) {
                        symmetry->for_each_distinct_bracelet_image(support,
                                [&](support::bitset image, size_t, bool, size_t) {
                            generator.add_forbidden_orbit(image);
                        });
                    } else {
                        generator.add_forbidden_orbit(support);
                    }
                } else if (result == support_result::candidate) {
                    if (symmetry) {
                        std::optional<candidate_type> solved_candidate;
                        if (conf_with_candidates_) solved_candidate = candidate_;

                        symmetry->for_each_distinct_bracelet_image(candidate_.support,
                                [&](support::bitset image, size_t multiplier_class, bool reflected, size_t right_shifts) {
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
                                    candidate_.vector[image_position] = solved_candidate->vector[position];
                                }
                            }
                            finalize_candidate(multiplier);
                        });
                    } else {
                        finalize_candidate(generator.add_forbidden_orbit(candidate_.support));
                    }
                }
                if (finishes_singleton_layer(support_size, 1)) {
                    probe_full_support_after_singletons(generator);
                    find_weakly_dominated_strategies();
                    if (!weakly_dominated_strategies.empty())
                        generator.add_forbidden_orbit(make_singleton(weakly_dominated_strategies.front()));
                }
            });

            if (conf_with_candidates_ && symmetry) {
                std::sort(candidates_.begin(), candidates_.end(), [](const candidate_type& left, const candidate_type& right) {
                    return left.support_size < right.support_size ||
                           (left.support_size == right.support_size && left.support < right.support);
                });
                for (size_t index = 0; index < candidates_.size(); ++index) candidates_[index].candidate_id = index + 1;
            }
        } else {
            support::NonCircularSupportGenerator generator(dimension_);
            // This lambda is the callback. The generator calls it once for each support and supplies both the mask and its size.
            generate_until_complete(generator, [&](support::bitset support, size_t support_size) {
                check_support_size_limit(support_size);
                if (support_size == dimension_) {
                    if (!full_support_checked) check_full_support_exact();
                    return;
                }
                log_support_size(support_size);
                const support_result result = analyze_support(support, support_size);
                if (result == support_result::curvature_ceiling) {
                    generator.add_forbidden(support);
                } else if (result == support_result::candidate) {
                    generator.add_forbidden(candidate_.support);
                    finalize_candidate(std::nullopt);
                }
                if (finishes_singleton_layer(support_size, dimension_)) {
                    probe_full_support_after_singletons(generator);
                    find_weakly_dominated_strategies();
                    for (const size_t strategy : weakly_dominated_strategies)
                        generator.add_forbidden(make_singleton(strategy));
                }
            });
        }
    } else {
        if (game_.compact_circular) {
            support::CircularSupportGeneratorMultiword generator(dimension_);
            std::optional<support::CircularAffineSymmetryMultiword> symmetry(std::in_place, game_.matrix);
            if (symmetry->multiplier_class_count() == 1) symmetry.reset();
            generate_until_complete(generator, [&](const support::bitset_multiword& support, size_t support_size) {
                check_support_size_limit(support_size);
                if (support_size == dimension_) {
                    if (!full_support_checked) check_full_support_exact();
                    return;
                }
                if (symmetry && !symmetry->is_representative(support)) return;
                log_support_size(support_size);
                const support_result result = analyze_support(support, support_size);
                if (result == support_result::curvature_ceiling) {
                    if (symmetry) {
                        symmetry->for_each_distinct_bracelet_image(support,
                                [&](const support::bitset_multiword& image, size_t, bool, size_t) {
                            generator.add_forbidden_orbit(image);
                        });
                    } else {
                        generator.add_forbidden_orbit(support);
                    }
                } else if (result == support_result::candidate) {
                    if (symmetry) {
                        std::optional<candidate_type> solved_candidate;
                        if (conf_with_candidates_) solved_candidate = candidate_;

                        symmetry->for_each_distinct_bracelet_image(candidate_.support,
                                [&](const support::bitset_multiword& image, size_t multiplier_class,
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
                                    candidate_.vector[image_position] = solved_candidate->vector[position];
                                }
                            }
                            finalize_candidate(multiplier);
                        });
                    } else {
                        finalize_candidate(generator.add_forbidden_orbit(candidate_.support));
                    }
                }
                if (finishes_singleton_layer(support_size, 1)) {
                    probe_full_support_after_singletons(generator);
                    find_weakly_dominated_strategies();
                    if (!weakly_dominated_strategies.empty())
                        generator.add_forbidden_orbit(make_singleton(weakly_dominated_strategies.front()));
                }
            });

            if (conf_with_candidates_ && symmetry) {
                std::sort(candidates_.begin(), candidates_.end(), [](const candidate_type& left, const candidate_type& right) {
                    return left.support_size < right.support_size ||
                           (left.support_size == right.support_size && left.support < right.support);
                });
                for (size_t index = 0; index < candidates_.size(); ++index) candidates_[index].candidate_id = index + 1;
            }
        } else {
            support::NonCircularSupportGeneratorMultiword generator(dimension_);
            generate_until_complete(generator, [&](const support::bitset_multiword& support, size_t support_size) {
                check_support_size_limit(support_size);
                if (support_size == dimension_) {
                    if (!full_support_checked) check_full_support_exact();
                    return;
                }
                log_support_size(support_size);
                const support_result result = analyze_support(support, support_size);
                if (result == support_result::curvature_ceiling) {
                    generator.add_forbidden(support);
                } else if (result == support_result::candidate) {
                    generator.add_forbidden(candidate_.support);
                    finalize_candidate(std::nullopt);
                }
                if (finishes_singleton_layer(support_size, dimension_)) {
                    probe_full_support_after_singletons(generator);
                    find_weakly_dominated_strategies();
                    for (const size_t strategy : weakly_dominated_strategies)
                        generator.add_forbidden(make_singleton(strategy));
                }
            });
        }
    }

    finish_log();
}

template<class SupportMask>
typename basic_analyzer<SupportMask>::support_result
basic_analyzer<SupportMask>::analyze_support(const SupportMask& support, size_t support_size)
{
    const bool candidate_found = exact_candidate_solver_.find(support, support_size, candidate_, conf_with_candidates_);
    if (!exact_candidate_solver_.reduced_hessian_is_negative_definite()) return support_result::curvature_ceiling;
    if (!candidate_found) return support_result::rejected;

    candidate_.support_size = support_size;
    candidate_.support = support;

    log_candidate();
    check_stability();
    return support_result::candidate;
}

template<class SupportMask>
void basic_analyzer<SupportMask>::finalize_candidate(std::optional<size_t> multiplier) {
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
