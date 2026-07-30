#include <mpfr.h>

#include <fracessa/candidate_rejector_double.hpp>

#include <algorithm>
#include <cfloat>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace candidate_rejection {
namespace {

constexpr size_t kMaxDimension = bs64::kMaxBitsetDimension;

constexpr bool has_binary64_format() noexcept
{
    return FLT_EVAL_METHOD == 0 &&
           sizeof(double) == sizeof(uint64_t) &&
           std::numeric_limits<double>::is_iec559 &&
           std::numeric_limits<double>::radix == 2 &&
           std::numeric_limits<double>::digits == 53 &&
           std::numeric_limits<double>::min_exponent == -1021 &&
           std::numeric_limits<double>::max_exponent == 1024 &&
           std::numeric_limits<double>::has_denorm == std::denorm_present;
}

// A strict nearest-rounded operation followed by one adjacent binary64 step
// encloses the exact real result. Non-finite endpoints invalidate the proof.
bool add_down(double left, double right, double& result) noexcept
{
    const double rounded = left + right;
    return round_down(rounded, result);
}

bool add_up(double left, double right, double& result) noexcept
{
    const double rounded = left + right;
    return round_up(rounded, result);
}

bool sub_down(double left, double right, double& result) noexcept
{
    const double rounded = left - right;
    return round_down(rounded, result);
}

bool sub_up(double left, double right, double& result) noexcept
{
    const double rounded = left - right;
    return round_up(rounded, result);
}

bool mul_down(double left, double right, double& result) noexcept
{
    const double rounded = left * right;
    return round_down(rounded, result);
}

bool mul_up(double left, double right, double& result) noexcept
{
    const double rounded = left * right;
    return round_up(rounded, result);
}

bool div_up(double numerator, double denominator, double& result) noexcept
{
    const double rounded = numerator / denominator;
    return round_up(rounded, result);
}

bool make_rational_enclosure(
    const fraction& value,
    mpfr_t converted,
    RationalEnclosure& result) noexcept
{
    // Both conversion stages use the same direction; converting one midpoint
    // and stepping it would not rigorously enclose an arbitrary rational.
    // FLINT reads this argument; the wrapper's escape hatch has no const overload.
    auto& raw_value = const_cast<fraction&>(value).data();
    fmpq_get_mpfr(converted, raw_value, MPFR_RNDD);
    const double lower = mpfr_get_d(converted, MPFR_RNDD);
    fmpq_get_mpfr(converted, raw_value, MPFR_RNDN);
    double midpoint = mpfr_get_d(converted, MPFR_RNDN);
    fmpq_get_mpfr(converted, raw_value, MPFR_RNDU);
    const double upper = mpfr_get_d(converted, MPFR_RNDU);

    if (!std::isfinite(lower) || !std::isfinite(midpoint) || !std::isfinite(upper) ||
        lower > midpoint || midpoint > upper) {
        return false;
    }

    double radius = 0.0;
    if (lower == upper) {
        midpoint = lower;
    } else {
        double lower_gap;
        double upper_gap;
        if (!sub_up(midpoint, lower, lower_gap) ||
            !sub_up(upper, midpoint, upper_gap)) {
            return false;
        }
        radius = std::max(lower_gap, upper_gap);
    }

    const double magnitude = std::max(std::abs(lower), std::abs(upper));
    if (!std::isfinite(radius) || radius < 0.0 || !std::isfinite(magnitude)) {
        return false;
    }

    result = {midpoint, radius, magnitude};
    return true;
}

} // namespace

bool triangular_solve_in_place(
    const linalg::matrix_dbl& matrix,
    size_t dimension,
    double* solution) noexcept
{
    // L has an implicit unit diagonal below the retained U triangle.
    for (size_t i = 0; i < dimension; ++i) {
        double sum = solution[i];
        for (size_t j = 0; j < i; ++j) {
            const double product = matrix(i, j) * solution[j];
            if (!std::isfinite(product)) return false;
            sum -= product;
            if (!std::isfinite(sum)) return false;
        }
        solution[i] = sum;
    }

    for (size_t i = dimension; i-- > 0;) {
        double sum = solution[i];
        for (size_t j = i + 1; j < dimension; ++j) {
            const double product = matrix(i, j) * solution[j];
            if (!std::isfinite(product)) return false;
            sum -= product;
            if (!std::isfinite(sum)) return false;
        }
        solution[i] = sum / matrix(i, i);
        if (!std::isfinite(solution[i])) return false;
    }
    return true;
}

namespace {

bool scan_outside_proposals(
    const linalg::matrix_dbl& game_midpoint,
    const uint8_t* support_indices,
    size_t support_size,
    const uint8_t* non_support_indices,
    size_t non_support_count,
    const double* solution,
    uint8_t* proposals,
    size_t& proposal_count) noexcept
{
    proposal_count = 0;
    const double payoff = solution[support_size];

    for (size_t i_pos = 0; i_pos < non_support_count; ++i_pos) {
        const size_t i = non_support_indices[i_pos];
        double rowsum = 0.0;
        for (size_t j_pos = 0; j_pos < support_size; ++j_pos) {
            const size_t j = support_indices[j_pos];
            const double product = game_midpoint(i, j) * solution[j_pos];
            if (!std::isfinite(product)) return false;
            rowsum += product;
            if (!std::isfinite(rowsum)) return false;
        }
        const double gain = rowsum - payoff;
        if (!std::isfinite(gain)) return false;
        if (gain > 0.0) proposals[proposal_count++] = static_cast<uint8_t>(i);
    }
    return true;
}

} // namespace

bool residual_enclosure(
    const linalg::matrix_dbl& game_midpoint,
    const linalg::matrix_dbl& game_radius,
    const uint8_t* support_indices,
    size_t support_size,
    const uint8_t* permutation,
    const double* solution,
    double* input_row_bound,
    double* residual_lower,
    double* residual_upper) noexcept
{
    // Enclose C0*z-b first, then widen by R_C*|z|. The final permutation
    // converts the original row order to the retained LU row order.
    const size_t dimension = support_size + 1;
    double unpermuted_lower[kMaxDimension];
    double unpermuted_upper[kMaxDimension];

    for (size_t row = 0; row < support_size; ++row) {
        const size_t i = support_indices[row];
        double lower = 0.0;
        double upper = 0.0;
        double input_radius = 0.0;
        double row_radius = 0.0;

        for (size_t column = 0; column < support_size; ++column) {
            const size_t j = support_indices[column];
            double product_lower;
            double product_upper;
            if (!mul_down(game_midpoint(i, j), solution[column], product_lower) ||
                !mul_up(game_midpoint(i, j), solution[column], product_upper) ||
                !add_down(lower, product_lower, lower) ||
                !add_up(upper, product_upper, upper)) {
                return false;
            }

            double radius_term;
            if (!mul_up(game_radius(i, j), std::abs(solution[column]), radius_term) ||
                !add_up(input_radius, radius_term, input_radius) ||
                !add_up(row_radius, game_radius(i, j), row_radius)) {
                return false;
            }
        }

        if (!sub_down(lower, solution[support_size], lower) ||
            !sub_up(upper, solution[support_size], upper) ||
            !sub_down(lower, input_radius, lower) ||
            !add_up(upper, input_radius, upper) ||
            lower > upper) {
            return false;
        }

        unpermuted_lower[row] = lower;
        unpermuted_upper[row] = upper;
        input_row_bound[row] = row_radius;
    }

    double lower = 0.0;
    double upper = 0.0;
    for (size_t column = 0; column < support_size; ++column) {
        if (!add_down(lower, solution[column], lower) ||
            !add_up(upper, solution[column], upper)) {
            return false;
        }
    }
    if (!sub_down(lower, 1.0, lower) ||
        !sub_up(upper, 1.0, upper) ||
        lower > upper) {
        return false;
    }
    unpermuted_lower[support_size] = lower;
    unpermuted_upper[support_size] = upper;
    input_row_bound[support_size] = 0.0;

    for (size_t row = 0; row < dimension; ++row) {
        residual_lower[row] = unpermuted_lower[permutation[row]];
        residual_upper[row] = unpermuted_upper[permutation[row]];
    }
    return true;
}

bool outside_gain_lower_bound(
    const linalg::matrix_dbl& game_midpoint,
    const linalg::matrix_dbl& game_radius,
    const linalg::matrix_dbl& game_magnitude,
    size_t outside_index,
    const uint8_t* support_indices,
    size_t support_size,
    const double* solution,
    double error,
    double& gain_lower) noexcept
{
    // |A*x-A0*x_hat| <= RA*|x_hat| + |A|*error; one more error covers payoff.
    double point_lower = 0.0;
    double gain_radius = 0.0;

    for (size_t j_pos = 0; j_pos < support_size; ++j_pos) {
        const size_t j = support_indices[j_pos];
        double product_lower;
        if (!mul_down(game_midpoint(outside_index, j), solution[j_pos], product_lower) ||
            !add_down(point_lower, product_lower, point_lower)) {
            return false;
        }

        double input_term;
        double solution_term;
        double term;
        if (!mul_up(game_radius(outside_index, j), std::abs(solution[j_pos]), input_term) ||
            !mul_up(game_magnitude(outside_index, j), error, solution_term) ||
            !add_up(input_term, solution_term, term) ||
            !add_up(gain_radius, term, gain_radius)) {
            return false;
        }
    }

    if (!sub_down(point_lower, solution[support_size], point_lower) ||
        !add_up(gain_radius, error, gain_radius)) {
        return false;
    }

    if (!sub_down(point_lower, gain_radius, gain_lower)) return false;
    return true;
}

bool round_down(double value, double& result) noexcept
{
    if (!std::isfinite(value)) return false;
    if (value == 0.0) {
        result = -std::numeric_limits<double>::denorm_min();
        return true;
    }

    uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    // Finite binary64 encodings are adjacent integers; negatives run backwards.
    if (value < 0.0) ++bits;
    else --bits;
    std::memcpy(&result, &bits, sizeof(result));
    return std::isfinite(result);
}

bool round_up(double value, double& result) noexcept
{
    if (!std::isfinite(value)) return false;
    if (value == 0.0) {
        result = std::numeric_limits<double>::denorm_min();
        return true;
    }

    uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    if (value > 0.0) ++bits;
    else --bits;
    std::memcpy(&result, &bits, sizeof(result));
    return std::isfinite(result);
}

const char* unavailable_reason() noexcept
{
#if FRACESSA_CANDIDATE_REJECTOR_DOUBLE_SUPPORTED
    if constexpr (!has_binary64_format()) {
        return "IEEE 754 binary64 arithmetic is unavailable";
    }
    if (std::fegetround() != FE_TONEAREST) {
        return "round-to-nearest is not active";
    }

    volatile double subnormal = std::numeric_limits<double>::denorm_min();
    volatile double two = 2.0;
    volatile double min_normal = std::numeric_limits<double>::min();
    volatile double half = 0.5;
    volatile double subnormal_input = subnormal * two;
    volatile double subnormal_output = min_normal * half;

    if (subnormal_input != 2.0 * std::numeric_limits<double>::denorm_min() ||
        subnormal_output != 0.5 * std::numeric_limits<double>::min()) {
        return "subnormal floating-point values are not preserved";
    }
    return nullptr;
#else
    return "this compiler is not supported by candidate-rejector-double";
#endif
}

bool rational_enclosure(const fraction& value, RationalEnclosure& result) noexcept
{
    mpfr_t converted;
    mpfr_init2(converted, 53);
    const bool success = make_rational_enclosure(value, converted, result);
    mpfr_clear(converted);
    return success;
}

bool factor_lu(
    linalg::matrix_dbl& matrix,
    size_t dimension,
    uint8_t* permutation) noexcept
{
    // Complete row swaps preserve earlier multipliers below U. permutation[i]
    // is the original row now occupying factor row i.
    for (size_t i = 0; i < dimension; ++i) {
        permutation[i] = static_cast<uint8_t>(i);
    }

    for (size_t k = 0; k < dimension; ++k) {
        size_t max_row = k;
        double max_value = std::abs(matrix(k, k));
        for (size_t i = k + 1; i < dimension; ++i) {
            const double value = std::abs(matrix(i, k));
            if (value > max_value) {
                max_value = value;
                max_row = i;
            }
        }
        if (max_value == 0.0 || !std::isfinite(max_value)) return false;

        if (max_row != k) {
            matrix.swap_rows(k, max_row);
            std::swap(permutation[k], permutation[max_row]);
        }

        const double pivot = matrix(k, k);
        for (size_t i = k + 1; i < dimension; ++i) {
            const double factor = matrix(i, k) / pivot;
            if (!std::isfinite(factor)) return false;
            matrix(i, k) = factor;
            for (size_t j = k + 1; j < dimension; ++j) {
                const double product = factor * matrix(k, j);
                if (!std::isfinite(product)) return false;
                const double updated = matrix(i, j) - product;
                if (!std::isfinite(updated)) return false;
                matrix(i, j) = updated;
            }
        }
    }
    return true;
}

bool absolute_lu_apply(
    const linalg::matrix_dbl& matrix,
    size_t dimension,
    const double* input,
    double* output) noexcept
{
    // Upward recurrences bound |U^-1|*|L^-1|*input without forming either inverse.
    double lower_solution[kMaxDimension];

    for (size_t i = 0; i < dimension; ++i) {
        if (!std::isfinite(input[i]) || input[i] < 0.0) return false;
        double sum = input[i];
        for (size_t j = 0; j < i; ++j) {
            double product;
            if (!mul_up(std::abs(matrix(i, j)), lower_solution[j], product) ||
                !add_up(sum, product, sum)) {
                return false;
            }
        }
        lower_solution[i] = sum;
    }

    for (size_t i = dimension; i-- > 0;) {
        double sum = lower_solution[i];
        for (size_t j = i + 1; j < dimension; ++j) {
            double product;
            if (!mul_up(std::abs(matrix(i, j)), output[j], product) ||
                !add_up(sum, product, sum)) {
                return false;
            }
        }
        const double pivot = std::abs(matrix(i, i));
        if (!(pivot > 0.0) || !std::isfinite(pivot) ||
            !div_up(sum, pivot, output[i])) {
            return false;
        }
    }
    return true;
}

bool prove_solution_error(
    const linalg::matrix_dbl& matrix,
    size_t dimension,
    const uint8_t* permutation,
    const double* input_row_bound,
    const double* residual_lower,
    const double* residual_upper,
    SolutionErrorBound& result) noexcept
{
    // Oishi-Rump's underflow-aware LU defect bound gives d >= |P*C-L*U|*1.
    // q < 1 then proves nonsingularity, and beta/(1-q) bounds solution error.
    constexpr double unit_roundoff = 0x1p-53;
    constexpr double underflow_unit = std::numeric_limits<double>::denorm_min();

    const double dimension_dbl = static_cast<double>(dimension);
    const double mu = dimension_dbl * unit_roundoff;
    const double rho = 1.0 - mu;
    double gamma;
    double delta;
    if (!(rho > 0.0) ||
        !div_up(mu, rho, gamma) ||
        !div_up(dimension_dbl, rho, delta)) {
        return false;
    }

    double upper_row_sum[kMaxDimension];
    for (size_t i = 0; i < dimension; ++i) {
        double row_sum = 0.0;
        for (size_t j = i; j < dimension; ++j) {
            if (!add_up(row_sum, std::abs(matrix(i, j)), row_sum)) return false;
        }
        upper_row_sum[i] = row_sum;
    }

    double defect[kMaxDimension];
    for (size_t i = 0; i < dimension; ++i) {
        double lu_row_sum = upper_row_sum[i];
        for (size_t j = 0; j < i; ++j) {
            double product;
            if (!mul_up(std::abs(matrix(i, j)), upper_row_sum[j], product) ||
                !add_up(lu_row_sum, product, lu_row_sum)) {
                return false;
            }
        }

        double factorization_term;
        double diagonal_term;
        double underflow_term;
        if (!mul_up(gamma, lu_row_sum, factorization_term) ||
            !add_up(dimension_dbl, std::abs(matrix(i, i)), diagonal_term) ||
            !mul_up(diagonal_term, underflow_unit, underflow_term) ||
            !mul_up(delta, underflow_term, underflow_term) ||
            !add_up(factorization_term, underflow_term, defect[i]) ||
            !add_up(defect[i], input_row_bound[permutation[i]], defect[i])) {
            return false;
        }
    }

    double propagated[kMaxDimension];
    if (!absolute_lu_apply(matrix, dimension, defect, propagated)) return false;
    double q = 0.0;
    for (size_t i = 0; i < dimension; ++i) q = std::max(q, propagated[i]);
    if (!std::isfinite(q) || !(q < 1.0)) return false;

    double residual_magnitude[kMaxDimension];
    for (size_t i = 0; i < dimension; ++i) {
        if (!std::isfinite(residual_lower[i]) || !std::isfinite(residual_upper[i]) ||
            residual_lower[i] > residual_upper[i]) {
            return false;
        }
        residual_magnitude[i] =
            std::max(std::abs(residual_lower[i]), std::abs(residual_upper[i]));
    }

    if (!absolute_lu_apply(matrix, dimension, residual_magnitude, propagated)) return false;
    double beta = 0.0;
    for (size_t i = 0; i < dimension; ++i) beta = std::max(beta, propagated[i]);

    double denominator;
    double error;
    if (!std::isfinite(beta) ||
        !sub_down(1.0, q, denominator) ||
        !(denominator > 0.0) ||
        !div_up(beta, denominator, error) ||
        error < 0.0) {
        return false;
    }

    result = {q, error};
    return true;
}

candidate_rejector_dbl::candidate_rejector_dbl(const linalg::matrix_frc& game_matrix)
    : game_frc_(game_matrix)
    , dimensions_(game_matrix.rows())
    , bounds_initialized_(false)
    , bounds_available_(false)
{
}

void candidate_rejector_dbl::initialize_bounds()
{
    bounds_initialized_ = true;
    bounds_available_ = false;

    if (unavailable_reason() || dimensions_ == 0) return;

    // Exact translation and positive scaling preserve every candidate inequality
    // while keeping all normalized matrix entries in [-1,1].
    const fraction translation = game_frc_(0, 0);
    fraction scale = fraction::zero();
    fraction difference;
    for (size_t i = 0; i < dimensions_; ++i) {
        for (size_t j = 0; j < dimensions_; ++j) {
            fraction::sub(difference, game_frc_(i, j), translation);
            if (difference.sgn() < 0) difference = -difference;
            if (difference > scale) scale = difference;
        }
    }
    if (scale.is_zero()) return;

    game_midpoint_ = linalg::matrix_dbl(dimensions_, dimensions_);
    game_radius_ = linalg::matrix_dbl(dimensions_, dimensions_);
    game_magnitude_ = linalg::matrix_dbl(dimensions_, dimensions_);

    fraction normalized;
    mpfr_t converted;
    mpfr_init2(converted, 53);
    for (size_t i = 0; i < dimensions_; ++i) {
        for (size_t j = 0; j < dimensions_; ++j) {
            fraction::sub(difference, game_frc_(i, j), translation);
            fraction::div(normalized, difference, scale);

            RationalEnclosure enclosure;
            if (!make_rational_enclosure(normalized, converted, enclosure)) {
                mpfr_clear(converted);
                return;
            }
            game_midpoint_(i, j) = enclosure.midpoint;
            game_radius_(i, j) = enclosure.radius;
            game_magnitude_(i, j) = enclosure.magnitude;
        }
    }
    mpfr_clear(converted);
    bounds_available_ = true;
}

bool candidate_rejector_dbl::proves_candidate_rejection(const bitset64& support, size_t support_size)
{
    if (!bounds_available_) {
        if (!bounds_initialized_) initialize_bounds();
        if (!bounds_available_) return false;
    }

    // Support, complement, solution, and proof vectors remain fixed stack storage.
    uint8_t support_indices[bs64::kMaxBitsetDimension];
    uint8_t non_support_indices[bs64::kMaxBitsetDimension];
    bs64::extract_set_indices(support, dimensions_, support_indices);
    const bitset64 complement = bs64::set_all_n_bits(dimensions_) & ~support;
    const size_t non_support_count =
        bs64::extract_set_indices(complement, dimensions_, non_support_indices);

    const size_t dimension = support_size + 1;
    if (compact_lu_.rows() != dimension) {
        compact_lu_ = linalg::matrix_dbl(dimension, dimension);
    }

    for (size_t row = 0; row < support_size; ++row) {
        const size_t i = support_indices[row];
        for (size_t column = 0; column < support_size; ++column) {
            compact_lu_(row, column) = game_midpoint_(i, support_indices[column]);
        }
        compact_lu_(row, support_size) = -1.0;
    }
    for (size_t column = 0; column < support_size; ++column) {
        compact_lu_(support_size, column) = 1.0;
    }
    compact_lu_(support_size, support_size) = 0.0;

    uint8_t permutation[bs64::kMaxBitsetDimension];
    if (!factor_lu(compact_lu_, dimension, permutation)) {
        return false;
    }

    double solution[bs64::kMaxBitsetDimension];
    for (size_t i = 0; i < dimension; ++i) {
        solution[i] = permutation[i] == support_size ? 1.0 : 0.0;
    }
    if (!triangular_solve_in_place(compact_lu_, dimension, solution)) {
        return false;
    }

    uint8_t support_proposals[bs64::kMaxBitsetDimension];
    size_t support_proposal_count = 0;
    for (size_t i = 0; i < support_size; ++i) {
        if (solution[i] <= 0.0) {
            support_proposals[support_proposal_count++] = static_cast<uint8_t>(i);
        }
    }

    uint8_t outside_proposals[bs64::kMaxBitsetDimension];
    size_t outside_proposal_count = 0;
    if (support_proposal_count == 0) {
        // Avoid the O(m^2) proof when the midpoint has no possible rejection.
        if (!scan_outside_proposals(
                game_midpoint_, support_indices, support_size,
                non_support_indices, non_support_count, solution,
                outside_proposals, outside_proposal_count)) {
            return false;
        }
        if (outside_proposal_count == 0) {
            return false;
        }
    }

    double input_row_bound[bs64::kMaxBitsetDimension];
    double residual_lower[bs64::kMaxBitsetDimension];
    double residual_upper[bs64::kMaxBitsetDimension];
    if (!residual_enclosure(
            game_midpoint_, game_radius_, support_indices, support_size,
            permutation, solution, input_row_bound,
            residual_lower, residual_upper)) {
        return false;
    }

    SolutionErrorBound proof;
    if (!prove_solution_error(
            compact_lu_, dimension, permutation, input_row_bound,
            residual_lower, residual_upper, proof)) {
        return false;
    }

    for (size_t i = 0; i < support_proposal_count; ++i) {
        double upper;
        if (!add_up(solution[support_proposals[i]], proof.error, upper)) {
            return false;
        }
        if (upper <= 0.0) {
            return true;
        }
    }

    if (support_proposal_count != 0) {
        if (!scan_outside_proposals(
                game_midpoint_, support_indices, support_size,
                non_support_indices, non_support_count, solution,
                outside_proposals, outside_proposal_count)) {
            return false;
        }
    }

    for (size_t i = 0; i < outside_proposal_count; ++i) {
        double gain_lower;
        if (!outside_gain_lower_bound(
                game_midpoint_, game_radius_, game_magnitude_, outside_proposals[i],
                support_indices, support_size, solution, proof.error, gain_lower)) {
            return false;
        }
        if (gain_lower > 0.0) {
            return true;
        }
    }

    return false;
}

} // namespace candidate_rejection
