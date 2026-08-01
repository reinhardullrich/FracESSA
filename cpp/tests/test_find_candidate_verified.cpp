#include <gmp.h>
#include <gtest/gtest.h>

#include <cfenv>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <fracessa/find_candidate_exact.hpp>
#include <fracessa/find_candidate_verified.hpp>
#include <fracessa/fracessa.hpp>

using namespace candidate_search;
using namespace linalg;

namespace {

fraction exact_double(double value)
{
    fraction result;
    mpq_t converted;
    mpq_init(converted);
    mpq_set_d(converted, value);
    fmpq_set_mpq(result.data(), converted);
    mpq_clear(converted);
    return result;
}

fraction absolute(fraction value)
{
    if (value.sgn() < 0) value = -value;
    return value;
}

fraction subtract(const fraction& left, const fraction& right)
{
    fraction result;
    fraction::sub(result, left, right);
    return result;
}

matrix_frc make_symmetric_2x2(const fraction& a, const fraction& b, const fraction& d)
{
    matrix_frc result(2, 2);
    result(0, 0) = a;
    result(0, 1) = b;
    result(1, 0) = b;
    result(1, 1) = d;
    return result;
}

class supported_find_candidate_verified_test : public ::testing::Test {
protected:
    void SetUp() override
    {
        if (unavailable_reason()) {
            GTEST_SKIP() << "Verified candidate search is unavailable";
        }
    }
};

} // namespace

TEST_F(supported_find_candidate_verified_test, StepsToAdjacentBinary64Values)
{
    double result;
    ASSERT_TRUE(round_down(1.0, result));
    EXPECT_EQ(result, std::nextafter(1.0, -std::numeric_limits<double>::infinity()));
    ASSERT_TRUE(round_up(1.0, result));
    EXPECT_EQ(result, std::nextafter(1.0, std::numeric_limits<double>::infinity()));
    ASSERT_TRUE(round_down(0.0, result));
    EXPECT_EQ(result, -std::numeric_limits<double>::denorm_min());
    ASSERT_TRUE(round_up(-0.0, result));
    EXPECT_EQ(result, std::numeric_limits<double>::denorm_min());
    EXPECT_FALSE(round_up(std::numeric_limits<double>::max(), result));
    EXPECT_FALSE(round_down(-std::numeric_limits<double>::max(), result));
}

TEST_F(supported_find_candidate_verified_test, DetectsUnsupportedRoundingMode)
{
    const int original_rounding = std::fegetround();
    ASSERT_EQ(std::fesetround(FE_UPWARD), 0);
    EXPECT_STREQ(unavailable_reason(), "round-to-nearest is not active");
    ASSERT_EQ(std::fesetround(original_rounding), 0);
}

TEST_F(supported_find_candidate_verified_test, RationalConversionEnclosesExactValues)
{
    rational_bounds enclosure;
    ASSERT_TRUE(get_rational_bounds(fraction(1, 2), enclosure));
    EXPECT_EQ(enclosure.midpoint, 0.5);
    EXPECT_EQ(enclosure.radius, 0.0);

    for (long numerator = -50; numerator <= 50; ++numerator) {
        for (long denominator = 1; denominator <= 47; ++denominator) {
            const fraction exact(
                static_cast<long long>(numerator),
                static_cast<long long>(denominator));
            ASSERT_TRUE(get_rational_bounds(exact, enclosure));

            const fraction midpoint = exact_double(enclosure.midpoint);
            const fraction radius = exact_double(enclosure.radius);
            EXPECT_LE(absolute(subtract(exact, midpoint)), radius);
            EXPECT_LE(absolute(exact), exact_double(enclosure.magnitude));
        }
    }

    ASSERT_TRUE(get_rational_bounds(fraction("1/1" + std::string(400, '0')), enclosure));
    EXPECT_GE(enclosure.radius, std::numeric_limits<double>::denorm_min());
    EXPECT_FALSE(get_rational_bounds(fraction("1" + std::string(400, '0')), enclosure));
}

TEST_F(supported_find_candidate_verified_test, RetainedLuSolvesAndTracksMultipleSwaps)
{
    matrix_dbl matrix(3, 3);
    matrix(0, 0) = 0.0; matrix(0, 1) = 1.0; matrix(0, 2) = 1.0;
    matrix(1, 0) = 2.0; matrix(1, 1) = 1.0; matrix(1, 2) = 0.0;
    matrix(2, 0) = 1.0; matrix(2, 1) = 3.0; matrix(2, 2) = 1.0;
    const matrix_dbl original = matrix;

    uint8_t permutation[bs64::kMaxBitsetDimension];
    ASSERT_TRUE(factor_lu(matrix, 3, permutation));
    EXPECT_EQ(permutation[0], 1);
    EXPECT_EQ(permutation[1], 2);
    EXPECT_EQ(permutation[2], 0);

    // The chosen exact solution is [2,3,1].
    const double right_hand_side[] = {4.0, 7.0, 12.0};
    double solution[bs64::kMaxBitsetDimension];
    for (size_t i = 0; i < 3; ++i) solution[i] = right_hand_side[permutation[i]];
    ASSERT_TRUE(triangular_solve_in_place(matrix, 3, solution));
    EXPECT_NEAR(solution[0], 2.0, 1e-14);
    EXPECT_NEAR(solution[1], 3.0, 1e-14);
    EXPECT_NEAR(solution[2], 1.0, 1e-14);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double reconstructed = 0.0;
            for (size_t k = 0; k < 3; ++k) {
                const double lower = i == k ? 1.0 : (i > k ? matrix(i, k) : 0.0);
                const double upper = k <= j ? matrix(k, j) : 0.0;
                reconstructed += lower * upper;
            }
            EXPECT_NEAR(reconstructed, original(permutation[i], j), 1e-14);
        }
    }
}

TEST_F(supported_find_candidate_verified_test, RetainedLuAllowsNegativeLastVariable)
{
    matrix_dbl matrix(3, 3);
    for (size_t i = 0; i < 3; ++i) matrix(i, i) = 1.0;

    uint8_t permutation[bs64::kMaxBitsetDimension];
    ASSERT_TRUE(factor_lu(matrix, 3, permutation));
    const double right_hand_side[] = {0.25, 0.75, -2.0};
    double solution[bs64::kMaxBitsetDimension];
    for (size_t i = 0; i < 3; ++i) solution[i] = right_hand_side[permutation[i]];
    ASSERT_TRUE(triangular_solve_in_place(matrix, 3, solution));
    EXPECT_EQ(solution[0], 0.25);
    EXPECT_EQ(solution[1], 0.75);
    EXPECT_EQ(solution[2], -2.0);
}

TEST_F(supported_find_candidate_verified_test, RetainedLuRejectsSingularAndNonFiniteWork)
{
    matrix_dbl singular(2, 2);
    singular(0, 0) = 1.0; singular(0, 1) = 1.0;
    singular(1, 0) = 2.0; singular(1, 1) = 2.0;
    uint8_t permutation[bs64::kMaxBitsetDimension];
    EXPECT_FALSE(factor_lu(singular, 2, permutation));

    matrix_dbl overflow(2, 2);
    overflow(0, 0) = std::numeric_limits<double>::denorm_min();
    overflow(0, 1) = std::numeric_limits<double>::max();
    overflow(1, 0) = std::numeric_limits<double>::denorm_min();
    overflow(1, 1) = -std::numeric_limits<double>::max();
    EXPECT_FALSE(factor_lu(overflow, 2, permutation));

    matrix_dbl non_finite(2, 2);
    non_finite(0, 0) = std::numeric_limits<double>::infinity();
    non_finite(1, 1) = 1.0;
    EXPECT_FALSE(factor_lu(non_finite, 2, permutation));
}

TEST_F(supported_find_candidate_verified_test, HandlesMaximumBorderedDimension)
{
    constexpr size_t dimension = bs64::kMaxBitsetDimension;
    matrix_dbl matrix(dimension, dimension);
    double right_hand_side[dimension];
    for (size_t i = 0; i < dimension; ++i) {
        matrix(i, i) = 1.0;
        right_hand_side[i] = static_cast<double>(i) - 31.0;
    }

    uint8_t permutation[dimension];
    ASSERT_TRUE(factor_lu(matrix, dimension, permutation));
    double solution[dimension];
    for (size_t i = 0; i < dimension; ++i) {
        solution[i] = right_hand_side[permutation[i]];
    }
    ASSERT_TRUE(triangular_solve_in_place(matrix, dimension, solution));
    for (size_t i = 0; i < dimension; ++i) EXPECT_EQ(solution[i], right_hand_side[i]);

    double zeros[dimension]{};
    solution_error_bound proof;
    ASSERT_TRUE(prove_solution_error(
        matrix, dimension, permutation, zeros, zeros, zeros, proof));
    EXPECT_LT(proof.q, 1.0);
    EXPECT_GE(proof.error, 0.0);
}

TEST_F(supported_find_candidate_verified_test, AbsoluteTriangularRecurrenceIsAnUpperBound)
{
    matrix_dbl matrix(3, 3);
    matrix(0, 0) = 2.0;  matrix(0, 1) = -1.0; matrix(0, 2) = 0.5;
    matrix(1, 0) = 0.25; matrix(1, 1) = 3.0; matrix(1, 2) = 1.0;
    matrix(2, 0) = -0.5; matrix(2, 1) = 0.75; matrix(2, 2) = 4.0;
    const double input[] = {0.5, 1.0, 2.0};
    double output[bs64::kMaxBitsetDimension];
    ASSERT_TRUE(absolute_lu_apply(matrix, 3, input, output));

    fraction y[3];
    for (size_t i = 0; i < 3; ++i) {
        y[i] = exact_double(input[i]);
        for (size_t j = 0; j < i; ++j) {
            y[i] += absolute(exact_double(matrix(i, j))) * y[j];
        }
    }
    fraction exact_output[3];
    for (size_t i = 3; i-- > 0;) {
        exact_output[i] = y[i];
        for (size_t j = i + 1; j < 3; ++j) {
            exact_output[i] += absolute(exact_double(matrix(i, j))) * exact_output[j];
        }
        exact_output[i] = exact_output[i] / absolute(exact_double(matrix(i, i)));
        EXPECT_FALSE(exact_double(output[i]) < exact_output[i]);
    }
}

TEST_F(supported_find_candidate_verified_test, ChoiceOneQBoundsExactTwoByTwoDefect)
{
    matrix_dbl original(2, 2);
    original(0, 0) = 3.0; original(0, 1) = 2.0;
    original(1, 0) = 7.0; original(1, 1) = 5.0;
    matrix_dbl matrix = original;
    uint8_t permutation[bs64::kMaxBitsetDimension];
    ASSERT_TRUE(factor_lu(matrix, 2, permutation));

    double zeros[bs64::kMaxBitsetDimension]{};
    solution_error_bound proof;
    ASSERT_TRUE(prove_solution_error(matrix, 2, permutation, zeros, zeros, zeros, proof));

    fraction product[2][2];
    product[0][0] = exact_double(matrix(0, 0));
    product[0][1] = exact_double(matrix(0, 1));
    product[1][0] = exact_double(matrix(1, 0)) * exact_double(matrix(0, 0));
    product[1][1] = exact_double(matrix(1, 0)) * exact_double(matrix(0, 1));
    product[1][1] += exact_double(matrix(1, 1));

    fraction defect[2][2];
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            defect[i][j] = subtract(
                exact_double(original(permutation[i], j)), product[i][j]);
        }
    }

    const fraction determinant = subtract(
        product[0][0] * product[1][1],
        product[0][1] * product[1][0]);
    fraction inverse[2][2];
    inverse[0][0] = product[1][1] / determinant;
    inverse[0][1] = -product[0][1] / determinant;
    inverse[1][0] = -product[1][0] / determinant;
    inverse[1][1] = product[0][0] / determinant;

    fraction exact_norm = fraction::zero();
    for (size_t i = 0; i < 2; ++i) {
        fraction row_sum = fraction::zero();
        for (size_t j = 0; j < 2; ++j) {
            fraction value = inverse[i][0] * defect[0][j];
            value += inverse[i][1] * defect[1][j];
            row_sum += absolute(value);
        }
        if (row_sum > exact_norm) exact_norm = row_sum;
    }
    EXPECT_FALSE(exact_double(proof.q) < exact_norm);
}

TEST_F(supported_find_candidate_verified_test, ResidualIntervalsContainExactResidual)
{
    matrix_dbl midpoint(3, 3);
    matrix_dbl radius(3, 3);
    midpoint(0, 0) = 1.0;  midpoint(0, 2) = -0.5;
    midpoint(2, 0) = 0.25; midpoint(2, 2) = 2.0;
    radius(0, 0) = 0.0625;   radius(0, 2) = 0.03125;
    radius(2, 0) = 0.015625; radius(2, 2) = 0.0625;

    const uint8_t support_indices[] = {0, 2};
    const uint8_t permutation[] = {2, 0, 1};
    const double solution[] = {0.25, 0.75, -0.5};
    double input_row_bound[bs64::kMaxBitsetDimension];
    double lower[bs64::kMaxBitsetDimension];
    double upper[bs64::kMaxBitsetDimension];
    ASSERT_TRUE(residual_enclosure(
        midpoint, radius, support_indices, 2, permutation, solution,
        input_row_bound, lower, upper));

    fraction exact_lower[3];
    fraction exact_upper[3];
    for (size_t row = 0; row < 2; ++row) {
        fraction residual = fraction::zero();
        fraction uncertainty = fraction::zero();
        for (size_t column = 0; column < 2; ++column) {
            residual += exact_double(midpoint(support_indices[row], support_indices[column])) *
                        exact_double(solution[column]);
            uncertainty += exact_double(radius(support_indices[row], support_indices[column])) *
                           absolute(exact_double(solution[column]));
        }
        residual -= exact_double(solution[2]);
        exact_lower[row] = residual;
        exact_lower[row] -= uncertainty;
        exact_upper[row] = residual;
        exact_upper[row] += uncertainty;

        fraction exact_row_radius = exact_double(radius(support_indices[row], support_indices[0]));
        exact_row_radius += exact_double(radius(support_indices[row], support_indices[1]));
        EXPECT_FALSE(exact_double(input_row_bound[row]) < exact_row_radius);
    }
    exact_lower[2] = exact_double(solution[0]);
    exact_lower[2] += exact_double(solution[1]);
    exact_lower[2] -= fraction::one();
    exact_upper[2] = exact_lower[2];

    for (size_t row = 0; row < 3; ++row) {
        const size_t original_row = permutation[row];
        EXPECT_LE(exact_double(lower[row]), exact_lower[original_row]);
        EXPECT_FALSE(exact_double(upper[row]) < exact_upper[original_row]);
    }
}

TEST_F(supported_find_candidate_verified_test, OutsideGainLowerBoundIsRigorous)
{
    matrix_dbl midpoint(2, 2);
    matrix_dbl radius(2, 2);
    matrix_dbl magnitude(2, 2);
    midpoint(1, 0) = 2.0;
    radius(1, 0) = 0.125;
    magnitude(1, 0) = 2.125;
    const uint8_t support_indices[] = {0};
    const double solution[] = {0.5, 0.25};
    const double error = 0.0625;

    double gain_lower;
    ASSERT_TRUE(outside_gain_lower_bound(
        midpoint, radius, magnitude, 1, support_indices, 1,
        solution, error, gain_lower));

    fraction exact_lower = subtract(
        exact_double(2.0) * exact_double(0.5), exact_double(0.25));
    exact_lower -= exact_double(0.125) * exact_double(0.5);
    exact_lower -= exact_double(2.125) * exact_double(error);
    exact_lower -= exact_double(error);
    EXPECT_LE(exact_double(gain_lower), exact_lower);
    EXPECT_GT(gain_lower, 0.0);
}

TEST_F(supported_find_candidate_verified_test, ProvesProfitableOutsideStrategy)
{
    const matrix_frc game = make_symmetric_2x2(fraction::zero(), fraction::one(), fraction::zero());
    find_candidate_verified verified(game);
    EXPECT_FALSE(verified.find(bitset64{1}, 1));
}

TEST_F(supported_find_candidate_verified_test, FallsBackAtOutsideGainBoundary)
{
    const matrix_frc game = make_symmetric_2x2(fraction::zero(), fraction::zero(), fraction::one());
    find_candidate_verified verified(game);
    EXPECT_TRUE(verified.find(bitset64{1}, 1));
}

TEST_F(supported_find_candidate_verified_test, ProvesNonPositiveSupportProbability)
{
    const matrix_frc game = make_symmetric_2x2(fraction::zero(), fraction::one(), fraction(3));
    find_candidate_verified verified(game);
    EXPECT_FALSE(verified.find(bitset64{3}, 2));
}

TEST_F(supported_find_candidate_verified_test, FallsBackAtZeroProbabilityBoundary)
{
    const matrix_frc game = make_symmetric_2x2(fraction::zero(), fraction::one(), fraction::one());
    find_candidate_verified verified(game);
    EXPECT_TRUE(verified.find(bitset64{3}, 2));
}

TEST_F(supported_find_candidate_verified_test, FallsBackForSingularAndNearSingularSystems)
{
    const matrix_frc singular_game = make_symmetric_2x2(fraction::zero(), fraction::one(), fraction::two());
    find_candidate_verified singular(singular_game);
    EXPECT_TRUE(singular.find(bitset64{3}, 2));

    const fraction epsilon("1/100000000000000000000");
    fraction nearly_two = fraction::two();
    nearly_two += epsilon;
    const matrix_frc nearly_singular_game = make_symmetric_2x2(fraction::zero(), fraction::one(), nearly_two);
    find_candidate_verified nearly_singular(nearly_singular_game);
    EXPECT_TRUE(nearly_singular.find(bitset64{3}, 2));
}

TEST_F(supported_find_candidate_verified_test, ReusesSameSizeScratchAfterFallback)
{
    matrix_frc game(3, 3);
    game(0, 0) = fraction::zero(); game(0, 1) = fraction::one(); game(0, 2) = fraction::one();
    game(1, 0) = fraction::one(); game(1, 1) = fraction::two(); game(1, 2) = fraction::zero();
    game(2, 0) = fraction::one(); game(2, 1) = fraction::zero(); game(2, 2) = fraction(3);
    find_candidate_verified verified(game);

    EXPECT_TRUE(verified.find(bitset64{3}, 2));
    EXPECT_FALSE(verified.find(bitset64{5}, 2));
}

TEST_F(supported_find_candidate_verified_test, NormalizationHandlesTinyScaleAndLargeTranslation)
{
    const fraction tiny("1/100000000000000000000");
    const matrix_frc tiny_game = make_symmetric_2x2(fraction::zero(), tiny, fraction::zero());
    find_candidate_verified tiny_scale(tiny_game);
    EXPECT_FALSE(tiny_scale.find(bitset64{1}, 1));

    const fraction base("100000000000000000000");
    fraction larger = base;
    larger += fraction::one();
    const matrix_frc translated_game = make_symmetric_2x2(base, larger, base);
    find_candidate_verified translated(translated_game);
    EXPECT_FALSE(translated.find(bitset64{1}, 1));
}

TEST(find_candidate_verified_fallback_test, ConstantMatrixFallsBack)
{
    const matrix_frc game = make_symmetric_2x2(fraction::one(), fraction::one(), fraction::one());
    find_candidate_verified constant(game);
    EXPECT_TRUE(constant.find(bitset64{1}, 1));
}

TEST(find_candidate_exact_test, BuildsVectorOnlyForRequestedSuccessfulCandidate)
{
    const matrix_frc game =
        make_symmetric_2x2(fraction::zero(), fraction::one(), fraction::zero());
    find_candidate_exact exact(game);

    candidate rejected;
    EXPECT_FALSE(exact.find(bitset64{1}, 1, rejected, true));
    EXPECT_EQ(rejected.vector.rows(), 0);

    candidate without_vector;
    ASSERT_TRUE(exact.find(bitset64{3}, 2, without_vector, false));
    EXPECT_EQ(without_vector.vector.rows(), 0);

    candidate with_vector;
    ASSERT_TRUE(exact.find(bitset64{3}, 2, with_vector, true));
    ASSERT_EQ(with_vector.vector.rows(), 2);
    ASSERT_EQ(with_vector.vector.cols(), 1);
    EXPECT_EQ(with_vector.vector(0, 0), fraction(1, 2));
    EXPECT_EQ(with_vector.vector(1, 0), fraction(1, 2));
}

TEST(find_candidate_verified_availability_test, UnavailableEnvironmentDisablesOnlyDefaultMode)
{
    const matrix_frc game =
        make_symmetric_2x2(fraction::zero(), fraction::one(), fraction::zero());
    const int original_rounding = std::fegetround();
    const bool runtime_failure_needed = unavailable_reason() == nullptr;
    if (runtime_failure_needed) ASSERT_EQ(std::fesetround(FE_UPWARD), 0);

    std::string message;
    try {
        ::fracessa analyzer(game, false);
    } catch (const std::exception& error) {
        message = error.what();
    }

    EXPECT_NO_THROW({
        ::fracessa analyzer(game, false, false, analysis_mode::exact);
    });
    EXPECT_NO_THROW({
        ::fracessa analyzer(game, false, false, analysis_mode::unsafe);
    });

    if (runtime_failure_needed) ASSERT_EQ(std::fesetround(original_rounding), 0);
    EXPECT_NE(message.find("--mode exact"), std::string::npos);
    EXPECT_NE(message.find("--mode unsafe"), std::string::npos);
}

TEST_F(supported_find_candidate_verified_test, NeverRejectsRandomExactCandidates)
{
    uint64_t state = 0x8f3d9b27a461ce55ULL;
    const auto next = [&state]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return state;
    };

    for (size_t sample = 0; sample < 120; ++sample) {
        const size_t dimension = 2 + sample % 5;
        matrix_frc game(dimension, dimension);
        for (size_t i = 0; i < dimension; ++i) {
            for (size_t j = i; j < dimension; ++j) {
                const long numerator = static_cast<long>(next() % 41) - 20;
                const long denominator = static_cast<long>(next() % 11) + 1;
                game(i, j) = game(j, i) = fraction(
                    static_cast<long long>(numerator),
                    static_cast<long long>(denominator));
            }
        }

        find_candidate_verified verified(game);
        find_candidate_exact exact(game);
        for (bitset64 support = 1; support < bs64::two_to_the_power_of(dimension); ++support) {
            const size_t support_size = bs64::count_set_bits(support);
            if (verified.find(support, support_size)) continue;

            candidate result;
            EXPECT_FALSE(exact.find(support, support_size, result, false))
                << "sample=" << sample << " support=" << support;
        }
    }
}
