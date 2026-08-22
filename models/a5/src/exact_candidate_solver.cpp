#include <fracessa/exact_candidate_solver.hpp>

#include <cassert>
#include <limits>
#include <stdexcept>

namespace fracessa::search {

namespace {

size_t checked_reduced_entry_cache_size(size_t dimension)
{
    const size_t maximum = std::numeric_limits<size_t>::max();
    if (dimension != 0 && dimension > maximum / dimension)
        throw std::overflow_error("Reduced-entry cache size overflow");
    const size_t square = dimension * dimension;
    if (dimension != 0 && square > maximum / dimension)
        throw std::overflow_error("Reduced-entry cache size overflow");
    return square * dimension;
}

} // namespace

exact_candidate_solver::exact_candidate_solver(const numeric::matrix_int& integer_game,
                                               const numeric::integer& game_denominator,
                                               support::SupportContext& support_context)
    : support_context_(support_context)
    , reduced_entry_cache_(checked_reduced_entry_cache_size(support_context.dimension()))
    , reduced_entry_cache_ready_(checked_reduced_entry_cache_size(support_context.dimension()), 0)
    , ffldlt_workspace_(support_context.dimension())
    , integer_game_(integer_game)
    , game_denominator_(game_denominator)
{
    support_indices_.reserve(support_context_.dimension());
    non_support_indices_.reserve(support_context_.dimension());
}

void exact_candidate_solver::resize_reduced_system(size_t reduced_dimension)
{
    if (reduced_system_dimension_ == reduced_dimension) return;

    reduced_system_dimension_ = reduced_dimension;
    reduced_system_.resize(reduced_system_dimension_, reduced_system_dimension_);
    right_hand_side_.resize(reduced_system_dimension_, 1);
    solution_numerators_.resize(reduced_system_dimension_, 1);
}

void exact_candidate_solver::build_reduced_system(const size_t* support_indices, size_t reduced_dimension)
{
    const size_t dimension = support_context_.dimension();
    const size_t reference = support_indices[0];
    const auto reference_diagonal = integer_game_(reference, reference);
    const size_t reference_cache_offset = reference * dimension * dimension;

    for (size_t row = 0; row < reduced_dimension; ++row) {
        const size_t i = support_indices[row + 1];
        auto row_offset = right_hand_side_(row, 0);
        row_offset.set_difference(reference_diagonal, integer_game_(i, reference));

        const size_t cache_row_offset = reference_cache_offset + i * dimension;
        auto* cached_values = reduced_entry_cache_.data() + cache_row_offset;
        auto* cached_ready = reduced_entry_cache_ready_.data() + cache_row_offset;

        for (size_t column = 0; column <= row; ++column) {
            const size_t j = support_indices[column + 1];
            auto& cached_value = cached_values[j];
            if (!cached_ready[j]) {
                cached_value.set_difference(integer_game_(i, j), integer_game_(reference, j));
                cached_value += row_offset;
                cached_ready[j] = 1;
            }
            reduced_system_(row, column) = cached_value;
        }
    }
}

numeric::integer::const_reference exact_candidate_solver::reduced_entry(size_t reference, size_t row, size_t column)
{
    const size_t dimension = support_context_.dimension();
    const size_t offset = reference * dimension * dimension + row * dimension + column;
    auto& value = reduced_entry_cache_[offset];
    if (!reduced_entry_cache_ready_[offset]) {
        value.set_difference(integer_game_(row, column), integer_game_(reference, column));
        value += integer_game_(reference, reference);
        value -= integer_game_(row, reference);
        reduced_entry_cache_ready_[offset] = 1;
    }
    return value;
}

size_t exact_candidate_solver::build_scaled_reduced_b(const support::Support& support,
                                                       const support::Support& outside_best_replies,
                                                       numeric::matrix_int& result)
{
    assert(!support_context_.empty(outside_best_replies));
    support_context_.extract_set_indices(support, support_indices_);
    support_context_.extract_set_indices(outside_best_replies, non_support_indices_);
    return build_scaled_reduced_b_from_indices(support_indices_.data(), support_indices_.size(),
                                               non_support_indices_.data(), non_support_indices_.size(), result);
}

size_t exact_candidate_solver::build_scaled_reduced_b_from_indices(const size_t* support_indices, size_t support_size,
                                                                   const size_t* outside_indices, size_t outside_size,
                                                                   numeric::matrix_int& result)
{
    assert(reduced_hessian_is_negative_definite_);
    assert(outside_size > 0);

    const size_t reduced_dimension = support_size - 1;
    const size_t reference = support_indices[0];

    if (reduced_dimension > 0) {
        assert(reduced_system_dimension_ == reduced_dimension);
        stability_solution_numerators_.resize(reduced_dimension, outside_size);

        for (size_t row = 0; row < reduced_dimension; ++row) {
            const size_t i = support_indices[row + 1];
            for (size_t column = 0; column < outside_size; ++column)
                stability_solution_numerators_(row, column) = reduced_entry(reference, i, outside_indices[column]);
        }

        ffldlt_workspace_.solve_factored_negative_definite_inplace(
            stability_solution_numerators_, solution_denominator_, reduced_system_);
    }

    result.resize(outside_size, outside_size);

    numeric::integer scaled_entry;
    for (size_t row = 0; row < outside_size; ++row) {
        for (size_t column = 0; column <= row; ++column) {
            scaled_entry.set_product(
                solution_denominator_, reduced_entry(reference, outside_indices[row], outside_indices[column]));

            for (size_t inner = 0; inner < reduced_dimension; ++inner) {
                const size_t i = support_indices[inner + 1];
                scaled_entry.submul(
                    reduced_entry(reference, i, outside_indices[row]), stability_solution_numerators_(inner, column));
            }

            scaled_entry.negate();
            result(row, column) = scaled_entry;
            if (row != column) result(column, row) = result(row, column);
        }
    }
    return outside_size;
}

void exact_candidate_solver::calculate_integer_payoff(numeric::integer& value, size_t strategy, size_t reference,
                                                       const size_t* support_indices, size_t reduced_dimension)
{
    value.set_product(integer_game_(strategy, reference), reference_numerator_);
    for (size_t position = 0; position < reduced_dimension; ++position)
        value.addmul(integer_game_(strategy, support_indices[position + 1]), solution_numerators_(position, 0));
}

void exact_candidate_solver::ensure_candidate_vector(candidate& result) const
{
    if (result.vector.size() != support_context_.dimension()) result.vector.resize(support_context_.dimension());
}

bool exact_candidate_solver::find(const support::Support& support, size_t support_size, candidate& result,
                                  bool materialize_output, support::Support& invaders)
{
    support_context_.extract_set_indices(support, support_indices_);
    support_context_.extract_unset_indices(support, non_support_indices_);
    return find_from_indices(support, support_size, result, materialize_output, support_indices_.data(), support_indices_.size(),
                             non_support_indices_.data(), non_support_indices_.size(), invaders);
}

bool exact_candidate_solver::find_from_indices(const support::Support& support, size_t support_size,
                                               candidate& result, bool materialize_output,
                                               const size_t* support_indices, size_t support_count,
                                               const size_t* non_support_indices, size_t non_support_count,
                                               support::Support& invaders)
{
    reduced_hessian_is_negative_definite_ = false;
    reduced_hessian_negative_inertia_ = 0;
    support_context_.clear(invaders);

    assert(support_count == support_size);
    static_cast<void>(support_size);

    const size_t reference = support_indices[0];
    const size_t reduced_dimension = support_count - 1;

    if (reduced_dimension == 0) {
        solution_denominator_.set_one();
        reference_numerator_.set_one();
        reduced_hessian_is_negative_definite_ = true;
    } else {
        resize_reduced_system(reduced_dimension);
        build_reduced_system(support_indices, reduced_dimension);

        numeric::fraction_free_ldlt_inertia inertia;
        const bool nonsingular = ffldlt_workspace_.solve_inplace(
            solution_numerators_, solution_denominator_, reduced_system_, right_hand_side_, inertia);
        reduced_hessian_negative_inertia_ = static_cast<size_t>(inertia.negative);
        if (!nonsingular) return false;
        reduced_hessian_is_negative_definite_ = reduced_hessian_negative_inertia_ == reduced_dimension;
        if (!reduced_hessian_is_negative_definite_) return false;

        reference_numerator_ = solution_denominator_;
        for (size_t position = 0; position < reduced_dimension; ++position) {
            const auto numerator = solution_numerators_(position, 0);
            if (numerator.sign() <= 0) return false;
            reference_numerator_ -= numerator;
        }
        if (reference_numerator_.sign() <= 0) return false;
    }

    calculate_integer_payoff(payoff_numerator_, reference, reference, support_indices, reduced_dimension);

    support_context_.copy(result.extended_support, support);
    bool has_invader = false;
    for (size_t position = 0; position < non_support_count; ++position) {
        const size_t outside_strategy = non_support_indices[position];
        calculate_integer_payoff(outside_payoff_numerator_, outside_strategy, reference, support_indices, reduced_dimension);
        const int comparison = outside_payoff_numerator_.compare(payoff_numerator_);
        if (comparison > 0) {
            has_invader = true;
            support_context_.set(invaders, outside_strategy);
        } else if (comparison == 0) {
            support_context_.set(result.extended_support, outside_strategy);
        }
    }
    if (has_invader) return false;
    result.extended_support_size = support_context_.count(result.extended_support);

    if (materialize_output) {
        payoff_denominator_.set_product(game_denominator_, solution_denominator_);
        result.payoff.set_ratio(payoff_numerator_, payoff_denominator_);
        result.payoff_dbl = result.payoff.to_dbl();

        ensure_candidate_vector(result);
        for (size_t position = 0; position < non_support_count; ++position)
            result.vector[non_support_indices[position]].set_zero();
        result.vector[reference].set_ratio(reference_numerator_, solution_denominator_);
        for (size_t position = 0; position < reduced_dimension; ++position)
            result.vector[support_indices[position + 1]].set_ratio(solution_numerators_(position, 0), solution_denominator_);
    }

    return true;
}

} // namespace fracessa::search
