#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

#include "ffldlt_candidate_solver.hpp"

#include <fracessa/find_candidate_exact.hpp>
#include <fracessa/matrix_parser.hpp>
#include <fracessa/supports.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;

const fmpq* raw_fraction(const fraction& value) noexcept
{
    return const_cast<fraction&>(value).data();
}

void ensure_candidate_vector(candidate& result, size_t dimension)
{
    if (result.vector.rows() != dimension || result.vector.cols() != 1) {
        result.vector = linalg::matrix_frc(dimension, 1);
    }
}

bool finish_rational_candidate(const linalg::matrix_frc& game, bitset64 support, size_t dimension,
                               const std::vector<fraction>& support_solution, candidate& result, bool materialize_vector)
{
    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    const size_t support_count = bs64::extract_set_indices(support, dimension, support_indices);
    const bitset64 complement = bs64::set_all_n_bits(dimension) & ~support;
    const size_t non_support_count = bs64::extract_set_indices(complement, dimension, non_support_indices);
    const size_t reference_index = support_indices[0];

    fraction payoff = fraction::zero();
    for (size_t position = 0; position < support_count; ++position) {
        payoff.addmul(game(reference_index, support_indices[position]), support_solution[position]);
    }

    result.extended_support = support;
    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        fraction rowsum = fraction::zero();
        for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
            rowsum.addmul(game(i, support_indices[j_pos]), support_solution[j_pos]);
        }
        if (rowsum > payoff) return false;
        if (rowsum == payoff) result.extended_support = bs64::set_bit_at_pos(result.extended_support, i);
    }

    result.payoff = payoff;
    result.payoff_dbl = payoff.to_dbl();
    result.extended_support_size = bs64::count_set_bits(result.extended_support);

    if (materialize_vector) {
        ensure_candidate_vector(result, dimension);
        for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
            result.vector(non_support_indices[i_pos], 0) = fraction::zero();
        }
        for (size_t i_pos = 0; i_pos < support_count; ++i_pos) {
            result.vector(support_indices[i_pos], 0) = support_solution[i_pos];
        }
    }
    return true;
}

/* The exact bordered Gaussian solver used immediately before the reduced-Hessian LDLT change. This is copied here so the historical algorithm
 * is built with the same compiler, dependencies, support generator, and timing harness as every new alternative. */
class old_bordered_gaussian {
public:
    explicit old_bordered_gaussian(const linalg::matrix_frc& game) noexcept
        : game_(game), dimension_(game.rows())
    {
    }

    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
    {
        const size_t system_dimension = support_size + 1;
        if (linear_system_.rows() != system_dimension) {
            linear_system_ = linalg::matrix_frc(system_dimension, system_dimension + 1);
        }

        uint8_t support_indices[bs64::kMaxBitsetDimension];
        uint8_t non_support_indices[bs64::kMaxBitsetDimension];
        const size_t support_count = bs64::extract_set_indices(support, dimension_, support_indices);
        const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
        const size_t non_support_count = bs64::extract_set_indices(complement, dimension_, non_support_indices);

        for (size_t row = 0; row < support_size; ++row) {
            const size_t i = support_indices[row];
            for (size_t column = 0; column < support_size; ++column) {
                linear_system_(row, column) = game_(i, support_indices[column]);
            }
            linear_system_(row, support_size) = fraction::neg_one();
            linear_system_(row, system_dimension) = fraction::zero();
        }
        for (size_t column = 0; column < support_size; ++column) {
            linear_system_(support_size, column) = fraction::one();
        }
        linear_system_(support_size, support_size) = fraction::zero();
        linear_system_(support_size, system_dimension) = fraction::one();

        auto& matrix = linear_system_;
        const size_t n = matrix.rows();
        for (size_t pivot_column = 0; pivot_column < n - 1; ++pivot_column) {
            size_t pivot_row = pivot_column;
            if (matrix(pivot_column, pivot_column).is_zero()) {
                bool found = false;
                for (size_t row = pivot_column + 1; row < n; ++row) {
                    if (!matrix(row, pivot_column).is_zero()) {
                        pivot_row = row;
                        found = true;
                        break;
                    }
                }
                if (!found) return false;
            }

            if (pivot_row != pivot_column) matrix.swap_rows(pivot_column, pivot_row);
            const fraction& pivot = matrix(pivot_column, pivot_column);

            for (size_t row = pivot_column + 1; row < n; ++row) {
                if (matrix(row, pivot_column).is_zero()) continue;

                fraction factor;
                fraction::div(factor, matrix(row, pivot_column), pivot);
                for (size_t column = pivot_column + 1; column <= n; ++column) {
                    matrix(row, column).submul(factor, matrix(pivot_column, column));
                }
                matrix(row, pivot_column) = fraction::zero();
            }
        }

        if (matrix(n - 1, n - 1).is_zero()) return false;

        linalg::matrix_frc solution(n, 1);
        for (size_t row = n; row-- > 0; ) {
            fraction sum = matrix(row, n);
            for (size_t column = row + 1; column < n; ++column) {
                sum.submul(matrix(row, column), solution(column, 0));
            }
            fraction::div(solution(row, 0), sum, matrix(row, row));
            if (row < n - 1 && solution(row, 0).sgn() <= 0) return false;
        }

        const fraction payoff = solution(support_size, 0);
        result.extended_support = support;
        for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
            const size_t i = non_support_indices[i_pos];
            fraction rowsum = fraction::zero();
            for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
                rowsum.addmul(game_(i, support_indices[j_pos]), solution(j_pos, 0));
            }
            if (rowsum > payoff) return false;
            if (rowsum == payoff) result.extended_support = bs64::set_bit_at_pos(result.extended_support, i);
        }

        result.payoff = payoff;
        result.payoff_dbl = payoff.to_dbl();
        result.extended_support_size = bs64::count_set_bits(result.extended_support);

        if (materialize_vector) {
            ensure_candidate_vector(result, dimension_);
            for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
                result.vector(non_support_indices[i_pos], 0) = fraction::zero();
            }
            for (size_t i_pos = 0; i_pos < support_count; ++i_pos) {
                result.vector(support_indices[i_pos], 0) = solution(i_pos, 0);
            }
        }
        return true;
    }

private:
    const linalg::matrix_frc& game_;
    size_t dimension_;
    linalg::matrix_frc linear_system_;
};

/*
 * Option A: clear each game row's denominators once, then solve the integer bordered system with FLINT's fraction-free LU.
 * If row i has positive denominator d_i and integer entries B_ij=d_i*A_ij, the candidate equations become
 *
 *     B_S*x - d_S*u = 0,       1^T*x = 1.
 *
 * FLINT returns one common solution denominator. Positivity and every outside-support comparison can therefore be decided
 * directly on the returned integer numerators; rational values are constructed only for a successful candidate.
 */
class integer_bordered_bareiss {
public:
    explicit integer_bordered_bareiss(const linalg::matrix_frc& game)
        : dimension_(game.rows())
    {
        fmpq_mat_t game_rational;
        fmpq_mat_init(game_rational, static_cast<slong>(dimension_), static_cast<slong>(dimension_));
        fmpz_mat_init(game_integer_, static_cast<slong>(dimension_), static_cast<slong>(dimension_));
        row_denominators_ = _fmpz_vec_init(static_cast<slong>(dimension_));

        for (size_t row = 0; row < dimension_; ++row) {
            for (size_t column = 0; column < dimension_; ++column) {
                fmpq_set(fmpq_mat_entry(game_rational, static_cast<slong>(row), static_cast<slong>(column)), raw_fraction(game(row, column)));
            }
        }
        fmpq_mat_get_fmpz_mat_rowwise(game_integer_, row_denominators_, game_rational);
        fmpq_mat_clear(game_rational);

        fmpz_mat_init(system_, 0, 0);
        fmpz_mat_init(right_hand_side_, 0, 0);
        fmpz_mat_init(solution_numerators_, 0, 0);
        fmpz_init(solution_denominator_);
        fmpz_init(outside_sum_);
        fmpz_init(outside_payoff_);
    }

    integer_bordered_bareiss(const integer_bordered_bareiss&) = delete;
    integer_bordered_bareiss& operator=(const integer_bordered_bareiss&) = delete;

    ~integer_bordered_bareiss()
    {
        fmpz_clear(outside_payoff_);
        fmpz_clear(outside_sum_);
        fmpz_clear(solution_denominator_);
        fmpz_mat_clear(solution_numerators_);
        fmpz_mat_clear(right_hand_side_);
        fmpz_mat_clear(system_);
        _fmpz_vec_clear(row_denominators_, static_cast<slong>(dimension_));
        fmpz_mat_clear(game_integer_);
    }

    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
    {
        uint8_t support_indices[bs64::kMaxBitsetDimension];
        uint8_t non_support_indices[bs64::kMaxBitsetDimension];
        const size_t support_count = bs64::extract_set_indices(support, dimension_, support_indices);
        const bitset64 complement = bs64::set_all_n_bits(dimension_) & ~support;
        const size_t non_support_count = bs64::extract_set_indices(complement, dimension_, non_support_indices);
        resize_system(support_size + 1);

        for (size_t row = 0; row < support_count; ++row) {
            const size_t i = support_indices[row];
            for (size_t column = 0; column < support_count; ++column) {
                const size_t j = support_indices[column];
                fmpz_set(fmpz_mat_entry(system_, static_cast<slong>(row), static_cast<slong>(column)),
                         fmpz_mat_entry(game_integer_, static_cast<slong>(i), static_cast<slong>(j)));
            }
            fmpz_neg(fmpz_mat_entry(system_, static_cast<slong>(row), static_cast<slong>(support_count)), row_denominators_ + i);
        }
        for (size_t column = 0; column < support_count; ++column) {
            fmpz_one(fmpz_mat_entry(system_, static_cast<slong>(support_count), static_cast<slong>(column)));
        }
        fmpz_zero(fmpz_mat_entry(system_, static_cast<slong>(support_count), static_cast<slong>(support_count)));

        if (!fmpz_mat_solve_fflu(solution_numerators_, solution_denominator_, system_, right_hand_side_)) return false;
        if (fmpz_sgn(solution_denominator_) < 0) {
            fmpz_mat_neg(solution_numerators_, solution_numerators_);
            fmpz_neg(solution_denominator_, solution_denominator_);
        }

        for (size_t position = 0; position < support_count; ++position) {
            if (fmpz_sgn(fmpz_mat_entry(solution_numerators_, static_cast<slong>(position), 0)) <= 0) return false;
        }

        const fmpz* payoff_numerator = fmpz_mat_entry(solution_numerators_, static_cast<slong>(support_count), 0);
        result.extended_support = support;
        for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
            const size_t i = non_support_indices[i_pos];
            fmpz_zero(outside_sum_);
            for (size_t j_pos = 0; j_pos < support_count; ++j_pos) {
                fmpz_addmul(outside_sum_, fmpz_mat_entry(game_integer_, static_cast<slong>(i), static_cast<slong>(support_indices[j_pos])),
                            fmpz_mat_entry(solution_numerators_, static_cast<slong>(j_pos), 0));
            }
            fmpz_mul(outside_payoff_, row_denominators_ + i, payoff_numerator);
            const int comparison = fmpz_cmp(outside_sum_, outside_payoff_);
            if (comparison > 0) return false;
            if (comparison == 0) result.extended_support = bs64::set_bit_at_pos(result.extended_support, i);
        }

        fmpq_set_fmpz_frac(result.payoff.data(), payoff_numerator, solution_denominator_);
        result.payoff_dbl = result.payoff.to_dbl();
        result.extended_support_size = bs64::count_set_bits(result.extended_support);

        if (materialize_vector) {
            ensure_candidate_vector(result, dimension_);
            for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
                result.vector(non_support_indices[i_pos], 0) = fraction::zero();
            }
            for (size_t position = 0; position < support_count; ++position) {
                fmpq_set_fmpz_frac(result.vector(support_indices[position], 0).data(),
                                   fmpz_mat_entry(solution_numerators_, static_cast<slong>(position), 0), solution_denominator_);
            }
        }
        return true;
    }

private:
    void resize_system(size_t system_dimension)
    {
        if (system_dimension_ == system_dimension) return;
        fmpz_mat_clear(solution_numerators_);
        fmpz_mat_clear(right_hand_side_);
        fmpz_mat_clear(system_);

        system_dimension_ = system_dimension;
        const slong size = static_cast<slong>(system_dimension_);
        fmpz_mat_init(system_, size, size);
        fmpz_mat_init(right_hand_side_, size, 1);
        fmpz_mat_init(solution_numerators_, size, 1);
        fmpz_mat_zero(right_hand_side_);
        fmpz_one(fmpz_mat_entry(right_hand_side_, size - 1, 0));
    }

    size_t dimension_;
    size_t system_dimension_ = 0;
    fmpz_mat_t game_integer_;
    fmpz* row_denominators_ = nullptr;
    fmpz_mat_t system_;
    fmpz_mat_t right_hand_side_;
    fmpz_mat_t solution_numerators_;
    fmpz_t solution_denominator_;
    fmpz_t outside_sum_;
    fmpz_t outside_payoff_;
};

/* Complete conservative hybrid: screen every support with integer FFLU, then call the unchanged current LDLT finder only for an exact candidate.
 * The second call repeats the candidate solve and outside checks, so this benchmark slightly overstates the integration cost. */
class integer_fflu_then_ldlt_candidates {
public:
    explicit integer_fflu_then_ldlt_candidates(const linalg::matrix_frc& game)
        : fflu_(game), ldlt_(game)
    {
    }

    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
    {
        if (!fflu_.find(support, support_size, screened_candidate_, false)) return false;
        if (!ldlt_.find(support, support_size, result, materialize_vector)) {
            throw std::runtime_error("integer FFLU accepted a candidate rejected by current LDLT");
        }
        return true;
    }

    bool reduced_hessian_is_negative_definite() const noexcept { return ldlt_.reduced_hessian_is_negative_definite(); }

private:
    integer_bordered_bareiss fflu_;
    candidate_search::find_candidate_exact ldlt_;
    candidate screened_candidate_;
};

/*
 * Option B: retain the reduced-Hessian equations, but hand the complete rational system to FLINT. FLINT clears its
 * denominators for every support and solves over the integers with fraction-free LU before returning rational y.
 */
class rational_fraction_free {
public:
    explicit rational_fraction_free(const linalg::matrix_frc& game)
        : game_(game), dimension_(game.rows()), support_solution_(dimension_)
    {
        fmpq_mat_init(reduced_hessian_, 0, 0);
        fmpq_mat_init(right_hand_side_, 0, 0);
        fmpq_mat_init(reduced_solution_, 0, 0);
    }

    rational_fraction_free(const rational_fraction_free&) = delete;
    rational_fraction_free& operator=(const rational_fraction_free&) = delete;

    ~rational_fraction_free()
    {
        fmpq_mat_clear(reduced_solution_);
        fmpq_mat_clear(right_hand_side_);
        fmpq_mat_clear(reduced_hessian_);
    }

    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector)
    {
        uint8_t support_indices[bs64::kMaxBitsetDimension];
        bs64::extract_set_indices(support, dimension_, support_indices);
        const size_t reduced_dimension = support_size - 1;

        if (reduced_dimension == 0) {
            support_solution_[0] = fraction::one();
        } else {
            resize_system(reduced_dimension);
            const size_t reference_index = support_indices[0];
            const fmpq* reference_diagonal = raw_fraction(game_(reference_index, reference_index));

            for (size_t row = 0; row < reduced_dimension; ++row) {
                const size_t i = support_indices[row + 1];
                fmpq* rhs = fmpq_mat_entry(right_hand_side_, static_cast<slong>(row), 0);
                fmpq_set(rhs, reference_diagonal);
                fmpq_sub(rhs, rhs, raw_fraction(game_(i, reference_index)));

                for (size_t column = 0; column <= row; ++column) {
                    const size_t j = support_indices[column + 1];
                    fmpq* value = fmpq_mat_entry(reduced_hessian_, static_cast<slong>(row), static_cast<slong>(column));
                    fmpq_set(value, raw_fraction(game_(i, j)));
                    fmpq_sub(value, value, raw_fraction(game_(i, reference_index)));
                    fmpq_sub(value, value, raw_fraction(game_(reference_index, j)));
                    fmpq_add(value, value, reference_diagonal);
                    if (row != column) {
                        fmpq_set(fmpq_mat_entry(reduced_hessian_, static_cast<slong>(column), static_cast<slong>(row)), value);
                    }
                }
            }

            if (!fmpq_mat_solve_fraction_free(reduced_solution_, reduced_hessian_, right_hand_side_)) return false;

            support_solution_[0] = fraction::one();
            for (size_t position = 0; position < reduced_dimension; ++position) {
                fmpq_set(support_solution_[position + 1].data(), fmpq_mat_entry(reduced_solution_, static_cast<slong>(position), 0));
                if (support_solution_[position + 1].sgn() <= 0) return false;
                support_solution_[0] -= support_solution_[position + 1];
            }
            if (support_solution_[0].sgn() <= 0) return false;
        }

        return finish_rational_candidate(game_, support, dimension_, support_solution_, result, materialize_vector);
    }

private:
    void resize_system(size_t reduced_dimension)
    {
        if (reduced_dimension_ == reduced_dimension) return;
        fmpq_mat_clear(reduced_solution_);
        fmpq_mat_clear(right_hand_side_);
        fmpq_mat_clear(reduced_hessian_);

        reduced_dimension_ = reduced_dimension;
        const slong size = static_cast<slong>(reduced_dimension_);
        fmpq_mat_init(reduced_hessian_, size, size);
        fmpq_mat_init(right_hand_side_, size, 1);
        fmpq_mat_init(reduced_solution_, size, 1);
    }

    const linalg::matrix_frc& game_;
    size_t dimension_;
    size_t reduced_dimension_ = 0;
    std::vector<fraction> support_solution_;
    fmpq_mat_t reduced_hessian_;
    fmpq_mat_t right_hand_side_;
    fmpq_mat_t reduced_solution_;
};

struct candidate_snapshot {
    bitset64 support = 0;
    bitset64 extended_support = 0;
    std::string payoff;
    std::vector<std::string> vector;
};

struct search_result {
    uint64_t visited_supports = 0;
    uint64_t candidate_count = 0;
    uint64_t checksum = 0;
    bool has_inertia = false;
    uint64_t negative_definite_count = 0;
    uint64_t inertia_checksum = 0;
    std::vector<candidate_snapshot> candidates;
    std::vector<uint8_t> negative_definite;
};

void mix_checksum(uint64_t& checksum, uint64_t value) noexcept
{
    checksum ^= value + 0x9e3779b97f4a7c15ULL + (checksum << 6) + (checksum >> 2);
}

candidate_snapshot make_snapshot(const candidate& value)
{
    candidate_snapshot snapshot;
    snapshot.support = value.support;
    snapshot.extended_support = value.extended_support;
    snapshot.payoff = value.payoff.to_string();
    snapshot.vector.reserve(value.vector.rows());
    for (size_t row = 0; row < value.vector.rows(); ++row) {
        snapshot.vector.push_back(value.vector(row, 0).to_string());
    }
    return snapshot;
}

template<bool provides_inertia, class Finder, class Generator>
search_result enumerate_candidates(Finder& finder, Generator& generator, bool capture_candidates)
{
    search_result search;
    search.has_inertia = provides_inertia;
    candidate current;
    generator.generate([&](bitset64 support, size_t support_size) {
        ++search.visited_supports;
        if (!finder.find(support, support_size, current, capture_candidates)) return;

        current.support = support;
        current.support_size = support_size;
        generator.add_forbidden(support);
        ++search.candidate_count;
        mix_checksum(search.checksum, support);
        mix_checksum(search.checksum, current.extended_support);
        if (capture_candidates) search.candidates.push_back(make_snapshot(current));

        if constexpr (provides_inertia) {
            const bool is_negative_definite = finder.reduced_hessian_is_negative_definite();
            search.negative_definite_count += is_negative_definite;
            mix_checksum(search.inertia_checksum, support);
            mix_checksum(search.inertia_checksum, is_negative_definite);
            if (capture_candidates) search.negative_definite.push_back(static_cast<uint8_t>(is_negative_definite));
        }
    });
    return search;
}

template<bool provides_inertia, class Finder>
search_result run_search(const linalg::matrix_frc& game, bool is_circular, bool capture_candidates)
{
    Finder finder(game);
    if (is_circular) {
        CircularSupportGenerator generator(game.rows());
        return enumerate_candidates<provides_inertia>(finder, generator, capture_candidates);
    }
    NonCircularSupportGenerator generator(game.rows());
    return enumerate_candidates<provides_inertia>(finder, generator, capture_candidates);
}

using runner = search_result (*)(const linalg::matrix_frc&, bool, bool);

search_result run_current(const linalg::matrix_frc& game, bool is_circular, bool capture_candidates)
{
    return run_search<true, candidate_search::find_candidate_exact>(game, is_circular, capture_candidates);
}

search_result run_old_gaussian(const linalg::matrix_frc& game, bool is_circular, bool capture_candidates)
{
    return run_search<false, old_bordered_gaussian>(game, is_circular, capture_candidates);
}

search_result run_integer_bareiss(const linalg::matrix_frc& game, bool is_circular, bool capture_candidates)
{
    return run_search<false, integer_bordered_bareiss>(game, is_circular, capture_candidates);
}

search_result run_integer_hybrid(const linalg::matrix_frc& game, bool is_circular, bool capture_candidates)
{
    return run_search<true, integer_fflu_then_ldlt_candidates>(game, is_circular, capture_candidates);
}

search_result run_ffldlt(const linalg::matrix_frc& game, bool is_circular, bool capture_candidates)
{
    return run_search<true, ffldlt_candidate_solver>(game, is_circular, capture_candidates);
}

search_result run_rational_fraction_free(const linalg::matrix_frc& game, bool is_circular, bool capture_candidates)
{
    return run_search<false, rational_fraction_free>(game, is_circular, capture_candidates);
}

bool same_snapshot(const candidate_snapshot& left, const candidate_snapshot& right)
{
    return left.support == right.support && left.extended_support == right.extended_support &&
           left.payoff == right.payoff && left.vector == right.vector;
}

void require_same_results(const search_result& expected, const search_result& actual, const std::string& label, bool compare_snapshots)
{
    if (expected.visited_supports != actual.visited_supports || expected.candidate_count != actual.candidate_count ||
        expected.checksum != actual.checksum) {
        throw std::runtime_error(label + " search result differs from current exact");
    }

    if (actual.has_inertia && (!expected.has_inertia || expected.negative_definite_count != actual.negative_definite_count ||
                               expected.inertia_checksum != actual.inertia_checksum)) {
        throw std::runtime_error(label + " reduced-Hessian inertia differs from current exact");
    }

    if (!compare_snapshots) return;
    if (expected.candidates.size() != actual.candidates.size()) {
        throw std::runtime_error(label + " candidate count differs from current exact");
    }
    for (size_t index = 0; index < expected.candidates.size(); ++index) {
        if (!same_snapshot(expected.candidates[index], actual.candidates[index])) {
            throw std::runtime_error(label + " candidate values differ at representative " + std::to_string(index));
        }
    }
    if (actual.has_inertia && expected.negative_definite != actual.negative_definite) {
        throw std::runtime_error(label + " per-candidate reduced-Hessian inertia differs from current exact");
    }
}

struct measurement {
    search_result result;
    uint64_t nanoseconds;
};

measurement measure(runner run, const linalg::matrix_frc& game, bool is_circular)
{
    const auto start = clock_type::now();
    search_result result = run(game, is_circular, false);
    const auto finish = clock_type::now();
    return {std::move(result), static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count())};
}

double median(std::vector<uint64_t> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    if (values.size() % 2 != 0) return static_cast<double>(values[middle]);
    return (static_cast<double>(values[middle - 1]) + static_cast<double>(values[middle])) / 2.0;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc < 2 || argc > 4) {
            std::cerr << "usage: exact_integer_solver_benchmark MATRIX [REPETITIONS | auto [TARGET_MS] | ffldlt_only [TARGET_MS]]\n";
            return 2;
        }

        const bool ffldlt_only = argc >= 3 && std::string(argv[2]) == "ffldlt_only";
        const bool automatic_sampling = argc >= 3 && std::string(argv[2]) == "auto";
        const bool adaptive_sampling = ffldlt_only || automatic_sampling;
        const size_t repetitions = adaptive_sampling ? 0 : (argc == 3 ? static_cast<size_t>(std::stoull(argv[2])) : 5);
        const uint64_t target_nanoseconds = adaptive_sampling
            ? (argc == 4 ? static_cast<uint64_t>(std::stoull(argv[3])) : 1000ULL) * 1'000'000ULL
            : 0;
        if (!adaptive_sampling && repetitions == 0) throw std::invalid_argument("repetitions must be positive");
        if (adaptive_sampling && target_nanoseconds == 0) throw std::invalid_argument("target milliseconds must be positive");

        linalg::matrix_frc game;
        bool is_circular = false;
        matrix_parser::parse_matrix_string(argv[1], game, is_circular);

        if (ffldlt_only) {
            std::vector<uint64_t> samples;
            const measurement first = measure(run_ffldlt, game, is_circular);
            samples.push_back(first.nanoseconds);
            const uint64_t first_nanoseconds = std::max<uint64_t>(first.nanoseconds, 1);
            const size_t sample_count = static_cast<size_t>(std::max<uint64_t>(1, (target_nanoseconds + first_nanoseconds - 1) / first_nanoseconds));
            for (size_t sample = 1; sample < sample_count; ++sample) {
                const measurement measured = measure(run_ffldlt, game, is_circular);
                require_same_results(first.result, measured.result, "ffldlt_candidate_solver", false);
                samples.push_back(measured.nanoseconds);
            }

            const double median_nanoseconds = median(samples);
            std::cout << "dimension\t" << game.rows() << '\n';
            std::cout << "circular\t" << (is_circular ? 1 : 0) << '\n';
            std::cout << "visited_supports\t" << first.result.visited_supports << '\n';
            std::cout << "candidate_representatives\t" << first.result.candidate_count << '\n';
            std::cout << "target_ms\t" << target_nanoseconds / 1'000'000ULL << '\n';
            std::cout << "variant\tsamples\tmedian_ns\tmedian_ms\n";
            std::cout << "ffldlt_candidate_solver\t" << samples.size() << '\t' << std::fixed << std::setprecision(0)
                      << median_nanoseconds << '\t' << std::setprecision(6) << median_nanoseconds / 1'000'000.0 << '\n';
            return 0;
        }

        constexpr size_t variant_count = 6;
        const std::array<std::string, variant_count> names{"current_reduced_ldlt", "old_bordered_gaussian", "integer_bordered_fflu",
                                                           "integer_fflu_then_ldlt_candidates", "ffldlt_candidate_solver", "rational_fraction_free"};
        const std::array<runner, variant_count> runners{run_current, run_old_gaussian, run_integer_bareiss, run_integer_hybrid,
                                                        run_ffldlt, run_rational_fraction_free};

        search_result expected;
        std::array<std::vector<uint64_t>, variant_count> samples;

        if (automatic_sampling) {
            // The detailed mode checks every candidate value. Avoid repeating multi-second cases in the full sweep:
            // use the first ordinary search as both the result check and first timing sample, then extend only short variants to the target time.
            for (size_t variant = 0; variant < runners.size(); ++variant) {
                const measurement first = measure(runners[variant], game, is_circular);
                if (variant == 0) expected = first.result;
                else require_same_results(expected, first.result, names[variant], false);
                samples[variant].push_back(first.nanoseconds);

                const uint64_t first_nanoseconds = std::max<uint64_t>(first.nanoseconds, 1);
                const size_t sample_count = static_cast<size_t>(std::max<uint64_t>(1, (target_nanoseconds + first_nanoseconds - 1) / first_nanoseconds));
                for (size_t sample = 1; sample < sample_count; ++sample) {
                    const measurement measured = measure(runners[variant], game, is_circular);
                    require_same_results(expected, measured.result, names[variant], false);
                    samples[variant].push_back(measured.nanoseconds);
                }
            }
        } else {
            expected = runners[0](game, is_circular, true);
            for (size_t variant = 1; variant < runners.size(); ++variant) {
                require_same_results(expected, runners[variant](game, is_circular, true), names[variant], true);
            }

            for (runner run : runners) require_same_results(expected, run(game, is_circular, false), "warm-up", false);

            for (size_t repetition = 0; repetition < repetitions; ++repetition) {
                for (size_t offset = 0; offset < runners.size(); ++offset) {
                    const size_t variant = (repetition + offset) % runners.size();
                    const measurement measured = measure(runners[variant], game, is_circular);
                    require_same_results(expected, measured.result, names[variant], false);
                    samples[variant].push_back(measured.nanoseconds);
                }
            }
        }

        std::array<double, variant_count> medians{};
        for (size_t variant = 0; variant < medians.size(); ++variant) medians[variant] = median(samples[variant]);

        std::cout << "dimension\t" << game.rows() << '\n';
        std::cout << "circular\t" << (is_circular ? 1 : 0) << '\n';
        std::cout << "visited_supports\t" << expected.visited_supports << '\n';
        std::cout << "candidate_representatives\t" << expected.candidate_count << '\n';
        if (automatic_sampling) std::cout << "target_ms\t" << target_nanoseconds / 1'000'000ULL << '\n';
        else std::cout << "repetitions\t" << repetitions << '\n';
        std::cout << "variant\tsamples\tmedian_ns\tmedian_ms\tchange_vs_old_gaussian_percent\tchange_vs_current_ldlt_percent\n";
        for (size_t variant = 0; variant < medians.size(); ++variant) {
            const double change_vs_old = 100.0 * (medians[variant] / medians[1] - 1.0);
            const double change_vs_current = 100.0 * (medians[variant] / medians[0] - 1.0);
            std::cout << names[variant] << '\t' << samples[variant].size() << '\t' << std::fixed << std::setprecision(0) << medians[variant] << '\t'
                      << std::setprecision(6) << medians[variant] / 1'000'000.0 << '\t' << std::setprecision(3) << change_vs_old << '\t'
                      << change_vs_current << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
