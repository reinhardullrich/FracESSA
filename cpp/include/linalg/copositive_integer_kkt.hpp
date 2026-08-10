#pragma once

/*
 * Isolated first-order KKT strict-copositivity checker adapted from Coposit's models/fracessa/solver.cpp. It is deliberately not
 * called by FracESSA's production stability path; the public function below exists for focused tests and later experiments.
 */

#include <fracessa/bitset64.hpp>
#include <fracessa/bitset_multiword.hpp>
#include <fracessa/support_generator_non_circular.hpp>
#include <linalg/copositive_integer.hpp>
#include <linalg/fraction_free_ldlt.hpp>
#include <linalg/matrix_integer.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace linalg {

namespace copositive_integer_kkt_detail {

// Private control-flow marker: an exact witness has settled the Boolean result. It never crosses the public function boundary.
struct not_strictly_copositive {};

inline bool contains(bitset64 support, size_t position) noexcept
{
    return bs64::is_set_at_pos(support, position);
}

inline bool contains(const bitset_multiword& support, size_t position) noexcept
{
    return support.is_set_at_pos(position);
}

template<class Index>
bool small_principal_matrix_is_strictly_copositive(const matrix_int& matrix, const Index* indices, size_t support_size)
{
    switch (support_size) {
    case 1:
        return is_strictly_copositive_1x1(matrix(indices[0], indices[0]));
    case 2:
        return is_strictly_copositive_2x2(
            matrix(indices[0], indices[0]), matrix(indices[0], indices[1]), matrix(indices[1], indices[1]));
    case 3:
        return is_strictly_copositive_3x3(
            matrix(indices[0], indices[0]), matrix(indices[0], indices[1]), matrix(indices[0], indices[2]),
            matrix(indices[1], indices[1]), matrix(indices[1], indices[2]), matrix(indices[2], indices[2]));
    default:
        assert(false && "small principal-matrix check requires support size one through three");
        return false;
    }
}

/*
 * Apply FracESSA's exact first-order candidate search to Q=-A. Every global maximum of x^T Q x is a KKT point. Therefore A is
 * strictly copositive exactly when every accepted KKT payoff is negative.
 */
class candidate_search {
public:
    explicit candidate_search(const matrix_int& matrix)
        : input_(matrix), game_(matrix), factorization_(matrix.rows())
    {
        game_.negate();
        support_indices_large_.reserve(matrix.rows());
    }

    bool solve()
    {
        try {
            if (game_.rows() <= bs64::kMaxBitsetDimension) solve_one_word();
            else solve_multiword();
        } catch (const not_strictly_copositive&) {
            return false;
        }

        if (!candidate_found_) throw std::runtime_error("KKT copositivity search found no KKT point");
        return true;
    }

private:
    void solve_one_word()
    {
        NonCircularSupportGenerator generator(game_.rows());
        generator.generate([&](bitset64 support, size_t support_size) {
            std::uint8_t support_indices[bs64::kMaxBitsetDimension];
            const size_t extracted_size = bs64::extract_set_indices(support, game_.rows(), support_indices);
            assert(extracted_size == support_size);
            if (check_support(support, support_indices, support_size)) generator.add_forbidden(support);
        });
    }

    void solve_multiword()
    {
        NonCircularSupportGeneratorMultiword generator(game_.rows());
        generator.generate([&](const bitset_multiword& support, size_t support_size) {
            support.extract_set_indices(support_indices_large_);
            assert(support_indices_large_.size() == support_size);
            if (check_support(support, support_indices_large_.data(), support_size)) generator.add_forbidden(support);
        });
    }

    template<class Support, class Index>
    bool check_support(const Support& support, const Index* support_indices, size_t support_size)
    {
        if (support_size <= 3 && !small_principal_matrix_is_strictly_copositive(input_, support_indices, support_size))
            throw not_strictly_copositive{};

        if (!find_candidate(support, support_indices, support_size)) return false;

        candidate_found_ = true;
        if (payoff_numerator_.sign() >= 0) throw not_strictly_copositive{};
        return true;
    }

    template<class Index>
    void build_reduced_system(const Index* support_indices, size_t support_size)
    {
        const size_t reduced_dimension = support_size - 1;
        reduced_system_.resize(reduced_dimension, reduced_dimension);
        solution_numerators_.resize(reduced_dimension, 1);

        const size_t reference = support_indices[0];
        const auto reference_diagonal = game_(reference, reference);
        for (size_t row = 0; row < reduced_dimension; ++row) {
            const size_t i = support_indices[row + 1];
            solution_numerators_(row, 0).set_difference(reference_diagonal, game_(i, reference));
            for (size_t column = 0; column <= row; ++column) {
                const size_t j = support_indices[column + 1];
                reduced_system_(row, column).set_difference(game_(i, j), game_(reference, j));
                reduced_system_(row, column) += solution_numerators_(row, 0);
            }
        }
    }

    template<class Index>
    void calculate_payoff(integer& payoff, size_t strategy, const Index* support_indices, size_t support_size)
    {
        payoff.set_product(game_(strategy, support_indices[0]), reference_numerator_);
        for (size_t position = 1; position < support_size; ++position)
            payoff.addmul(game_(strategy, support_indices[position]), solution_numerators_(position - 1, 0));
    }

    template<class Support, class Index>
    bool find_candidate(const Support& support, const Index* support_indices, size_t support_size)
    {
        if (support_size == 1) {
            solution_denominator_.set_one();
            reference_numerator_.set_one();
        } else {
            build_reduced_system(support_indices, support_size);
            if (!factorization_.factorize_inplace(reduced_system_)) return false;
            factorization_.solve_inplace(solution_numerators_, solution_denominator_, reduced_system_);
            assert(solution_denominator_.sign() > 0);

            reference_numerator_ = solution_denominator_;
            for (size_t position = 1; position < support_size; ++position) {
                const auto numerator = solution_numerators_(position - 1, 0);
                if (numerator.sign() <= 0) return false;
                reference_numerator_ -= numerator;
            }
            if (reference_numerator_.sign() <= 0) return false;
        }

        calculate_payoff(payoff_numerator_, support_indices[0], support_indices, support_size);
        for (size_t strategy = 0; strategy < game_.rows(); ++strategy) {
            if (contains(support, strategy)) continue;
            calculate_payoff(outside_payoff_numerator_, strategy, support_indices, support_size);
            if (outside_payoff_numerator_.compare(payoff_numerator_) > 0) return false;
        }
        return true;
    }

    const matrix_int& input_;
    matrix_int game_;
    fraction_free_ldlt_factorization factorization_;
    matrix_int reduced_system_;
    matrix_int solution_numerators_;
    integer solution_denominator_;
    integer reference_numerator_;
    integer payoff_numerator_;
    integer outside_payoff_numerator_;
    std::vector<size_t> support_indices_large_;
    bool candidate_found_ = false;
};

} // namespace copositive_integer_kkt_detail

/*
 * Experimental exact checker. It is callable directly but is not connected to FracESSA's runtime copositivity dispatch.
 */
inline bool is_strictly_copositive_kkt(const matrix_int& matrix)
{
    const size_t dimension = matrix.rows();
    if (dimension == 0 || matrix.cols() != dimension) throw std::invalid_argument("matrix must be nonempty and square");
    for (size_t row = 0; row < dimension; ++row) {
        for (size_t column = row + 1; column < dimension; ++column) {
            if (matrix(row, column).compare(matrix(column, row)) != 0)
                throw std::invalid_argument("matrix must be symmetric");
        }
    }

    if (dimension == 1) return is_strictly_copositive_1x1(matrix(0, 0));
    if (dimension == 2) return is_strictly_copositive_2x2(matrix(0, 0), matrix(0, 1), matrix(1, 1));
    if (dimension == 3)
        return is_strictly_copositive_3x3(
            matrix(0, 0), matrix(0, 1), matrix(0, 2), matrix(1, 1), matrix(1, 2), matrix(2, 2));

    return copositive_integer_kkt_detail::candidate_search(matrix).solve();
}

} // namespace linalg
