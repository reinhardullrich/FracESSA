#include <gtest/gtest.h>

#include <cstdint>

#include <fracessa/find_candidate_safe.hpp>
#include <linalg/fraction_free_ldlt.hpp>

/*
 * Exact reduced-Hessian LDL^T tests through the candidate procedure that owns it.
 */

TEST(LinearSolverFractionTest, SimpleCandidate) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::zero(); game(0, 1) = fraction::one();
    game(1, 0) = fraction::one();  game(1, 1) = fraction::zero();

    candidate result;
    candidate_search::find_candidate_safe finder(game);
    ASSERT_TRUE(finder.find(0b11, 2, result, true));
    EXPECT_EQ(result.vector(0, 0), fraction(1, 2));
    EXPECT_EQ(result.vector(1, 0), fraction(1, 2));
    EXPECT_EQ(result.payoff, fraction(1, 2));
    EXPECT_TRUE(finder.reduced_hessian_is_negative_definite());
}

TEST(LinearSolverFractionTest, AllowsNegativePayoffVariable) {
    linalg::matrix_frc game(1, 1);
    game(0, 0) = fraction(-2);

    candidate result;
    candidate_search::find_candidate_safe finder(game);
    ASSERT_TRUE(finder.find(0b1, 1, result, false));
    EXPECT_EQ(result.payoff, fraction(-2));
    EXPECT_TRUE(finder.reduced_hessian_is_negative_definite());
}

TEST(LinearSolverFractionTest, HandlesNonsingularZeroDiagonalTwoByTwoPivot) {
    /*
     * With strategy 0 as the reference, this game reduces exactly to
     *
     *     H = [0 1],      r = [1/4].
     *         [1 0]           [1/4]
     *
     * H is nonsingular but has no nonzero diagonal pivot. A scalar-only LDL^T would incorrectly call it singular; the exact 2x2 pivot must instead
     * recover y=(1/4,1/4), hence x=(1/2,1/4,1/4).
     */
    linalg::matrix_frc game(3, 3);
    game(0, 0) = fraction::zero(); game(0, 1) = fraction(-1, 4); game(0, 2) = fraction(-1, 4);
    game(1, 0) = fraction(-1, 4);  game(1, 1) = fraction(-1, 2); game(1, 2) = fraction(1, 2);
    game(2, 0) = fraction(-1, 4);  game(2, 1) = fraction(1, 2);  game(2, 2) = fraction(-1, 2);

    candidate result;
    candidate_search::find_candidate_safe finder(game);
    ASSERT_TRUE(finder.find(0b111, 3, result, true));
    EXPECT_EQ(result.vector(0, 0), fraction(1, 2));
    EXPECT_EQ(result.vector(1, 0), fraction(1, 4));
    EXPECT_EQ(result.vector(2, 0), fraction(1, 4));
    EXPECT_EQ(result.payoff, fraction(-1, 8));

    // The 2x2 pivot has one eigenvalue of each sign, so H is indefinite.
    EXPECT_FALSE(finder.reduced_hessian_is_negative_definite());
}

TEST(LinearSolverFractionTest, RejectsNonPositiveSupportVariable) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::one(); game(0, 1) = fraction::zero();
    game(1, 0) = fraction::zero(); game(1, 1) = fraction::zero();

    candidate result;
    candidate_search::find_candidate_safe finder(game);
    EXPECT_FALSE(finder.find(0b11, 2, result, false));
    EXPECT_FALSE(finder.reduced_hessian_is_negative_definite());
}

TEST(LinearSolverFractionTest, RejectsSingularSystem) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::one(); game(0, 1) = fraction::one();
    game(1, 0) = fraction::one(); game(1, 1) = fraction::one();

    candidate result;
    candidate_search::find_candidate_safe finder(game);
    EXPECT_FALSE(finder.find(0b11, 2, result, false));
}

TEST(LinearSolverFractionTest, ReusesNegativeDefiniteFactorForMultipleRightHandSides) {
    linalg::matrix_int system(3, 3);
    linalg::integer value(-4);
    system(0, 0) = value;
    value = linalg::integer(-1);
    system(1, 0) = value;
    value = linalg::integer(-3);
    system(1, 1) = value;
    value = linalg::integer(-1);
    system(2, 0) = value;
    system(2, 1) = value;
    value = linalg::integer(-2);
    system(2, 2) = value;

    linalg::matrix_int initial_right_hand_side(3, 1);
    linalg::matrix_int initial_solution(3, 1);
    linalg::integer denominator;
    linalg::fraction_free_ldlt_inertia inertia;
    linalg::kkt_fraction_free_ldlt_workspace workspace(3);

    ASSERT_EQ(workspace.solve_inplace(initial_solution, denominator, system, initial_right_hand_side, inertia), 1);
    EXPECT_EQ(inertia.positive, 0);
    EXPECT_EQ(inertia.negative, 3);

    // H*X for X columns (1,3,5) and (2,4,6). The exact solution is integral, while the solver deliberately retains |det(H)|=17.
    linalg::matrix_int right_hand_sides(3, 2);
    value = linalg::integer(-12);
    right_hand_sides(0, 0) = value;
    value = linalg::integer(-18);
    right_hand_sides(0, 1) = value;
    value = linalg::integer(-15);
    right_hand_sides(1, 0) = value;
    value = linalg::integer(-20);
    right_hand_sides(1, 1) = value;
    value = linalg::integer(-14);
    right_hand_sides(2, 0) = value;
    value = linalg::integer(-18);
    right_hand_sides(2, 1) = value;

    workspace.solve_factored_negative_definite_inplace(right_hand_sides, denominator, system);

    value = linalg::integer(17);
    EXPECT_EQ(denominator.compare(value), 0);
    value = linalg::integer(17);
    EXPECT_EQ(right_hand_sides(0, 0).compare(value), 0);
    value = linalg::integer(34);
    EXPECT_EQ(right_hand_sides(0, 1).compare(value), 0);
    value = linalg::integer(51);
    EXPECT_EQ(right_hand_sides(1, 0).compare(value), 0);
    value = linalg::integer(68);
    EXPECT_EQ(right_hand_sides(1, 1).compare(value), 0);
    value = linalg::integer(85);
    EXPECT_EQ(right_hand_sides(2, 0).compare(value), 0);
    value = linalg::integer(102);
    EXPECT_EQ(right_hand_sides(2, 1).compare(value), 0);
}

TEST(FractionFreeLDLTFactorizationTest, RetainsFactorizationForSeveralRightHandSides) {
    linalg::matrix_int system(2, 2);
    system(0, 0) = linalg::integer(4);
    system(1, 0) = linalg::integer(1);
    system(1, 1) = linalg::integer(3);

    linalg::fraction_free_ldlt_factorization factorization(2);
    ASSERT_EQ(factorization.factorize_inplace(system), 1);
    EXPECT_EQ(factorization.determinant().compare(linalg::integer(11)), 0);
    EXPECT_TRUE(factorization.is_positive_definite());

    // The columns are b=(1,2) and b=(7,-1), whose solutions are (1/11,7/11) and (2,-1).
    linalg::matrix_int right_hand_sides(2, 2);
    right_hand_sides(0, 0) = linalg::integer(1);
    right_hand_sides(1, 0) = linalg::integer(2);
    right_hand_sides(0, 1) = linalg::integer(7);
    right_hand_sides(1, 1) = linalg::integer(-1);
    linalg::integer denominator;

    factorization.solve_inplace(right_hand_sides, denominator, system);

    EXPECT_EQ(denominator.compare(linalg::integer(11)), 0);
    EXPECT_EQ(right_hand_sides(0, 0).compare(linalg::integer(1)), 0);
    EXPECT_EQ(right_hand_sides(1, 0).compare(linalg::integer(7)), 0);
    EXPECT_EQ(right_hand_sides(0, 1).compare(linalg::integer(22)), 0);
    EXPECT_EQ(right_hand_sides(1, 1).compare(linalg::integer(-11)), 0);
}

TEST(FractionFreeLDLTFactorizationTest, ReplaysZeroDiagonalCoordinateOperations) {
    // After eliminating coordinate 0, the remaining 2x2 block has zero diagonal and a nonzero off-diagonal entry.
    linalg::matrix_int system(3, 3);
    system(0, 0) = linalg::integer(1);
    system(1, 0) = linalg::integer(1);
    system(1, 1) = linalg::integer(1);
    system(2, 0) = linalg::integer(0);
    system(2, 1) = linalg::integer(1);
    system(2, 2) = linalg::integer(0);

    linalg::fraction_free_ldlt_factorization factorization(3);
    ASSERT_EQ(factorization.factorize_inplace(system), 1);
    EXPECT_EQ(factorization.determinant().compare(linalg::integer(-1)), 0);
    EXPECT_FALSE(factorization.is_positive_definite());

    // The columns are b=(1,2,3) and b=(1,5,-1), whose solutions are (-2,3,1) and (2,-1,4).
    linalg::matrix_int right_hand_sides(3, 2);
    right_hand_sides(0, 0) = linalg::integer(1);
    right_hand_sides(1, 0) = linalg::integer(2);
    right_hand_sides(2, 0) = linalg::integer(3);
    right_hand_sides(0, 1) = linalg::integer(1);
    right_hand_sides(1, 1) = linalg::integer(5);
    right_hand_sides(2, 1) = linalg::integer(-1);
    linalg::integer denominator;

    factorization.solve_inplace(right_hand_sides, denominator, system);

    EXPECT_EQ(denominator.compare(linalg::integer(1)), 0);
    EXPECT_EQ(right_hand_sides(0, 0).compare(linalg::integer(-2)), 0);
    EXPECT_EQ(right_hand_sides(1, 0).compare(linalg::integer(3)), 0);
    EXPECT_EQ(right_hand_sides(2, 0).compare(linalg::integer(1)), 0);
    EXPECT_EQ(right_hand_sides(0, 1).compare(linalg::integer(2)), 0);
    EXPECT_EQ(right_hand_sides(1, 1).compare(linalg::integer(-1)), 0);
    EXPECT_EQ(right_hand_sides(2, 1).compare(linalg::integer(4)), 0);
}

TEST(FractionFreeLDLTFactorizationTest, ReturnsZeroDeterminantForSingularMatrix) {
    linalg::matrix_int system(2, 2);
    system(0, 0) = linalg::integer(1);
    system(1, 0) = linalg::integer(1);
    system(1, 1) = linalg::integer(1);

    linalg::fraction_free_ldlt_factorization factorization(2);
    EXPECT_EQ(factorization.factorize_inplace(system), 0);
    EXPECT_TRUE(factorization.determinant().is_zero());
    EXPECT_FALSE(factorization.is_positive_definite());
}

TEST(FractionFreeLDLTFactorizationTest, SolvesWithArbitraryPrecisionEntries) {
    linalg::integer big;
    ASSERT_EQ(fmpz_set_str(big.native_handle(), "123456789012345678901234567890", 10), 0);
    const auto set_multiple = [&big](linalg::integer::reference destination, ulong multiplier) {
        fmpz_mul_ui(destination.native_handle(), big.native_handle(), multiplier);
    };

    // A=big*[[4,1,2],[1,3,1],[2,1,5]] forces arbitrary-precision multiply/subtract and exact division steps.
    linalg::matrix_int system(3, 3);
    set_multiple(system(0, 0), 4);
    set_multiple(system(1, 0), 1);
    set_multiple(system(1, 1), 3);
    set_multiple(system(2, 0), 2);
    set_multiple(system(2, 1), 1);
    set_multiple(system(2, 2), 5);

    linalg::integer expected_determinant;
    fmpz_pow_ui(expected_determinant.native_handle(), big.native_handle(), 3);
    fmpz_mul_ui(expected_determinant.native_handle(), expected_determinant.native_handle(), 43);

    linalg::fraction_free_ldlt_factorization factorization(3);
    ASSERT_EQ(factorization.factorize_inplace(system), 1);
    EXPECT_EQ(factorization.determinant().compare(expected_determinant), 0);
    EXPECT_TRUE(factorization.is_positive_definite());

    // b=A*(1,2,3)=big*(12,10,19).
    linalg::matrix_int right_hand_side(3, 1);
    set_multiple(right_hand_side(0, 0), 12);
    set_multiple(right_hand_side(1, 0), 10);
    set_multiple(right_hand_side(2, 0), 19);
    linalg::integer denominator;

    factorization.solve_inplace(right_hand_side, denominator, system);

    EXPECT_EQ(denominator.compare(expected_determinant), 0);
    EXPECT_EQ(right_hand_side(0, 0).compare(expected_determinant), 0);
    linalg::integer twice_determinant(expected_determinant);
    twice_determinant.multiply(2);
    EXPECT_EQ(right_hand_side(1, 0).compare(twice_determinant), 0);
    linalg::integer three_times_determinant(expected_determinant);
    three_times_determinant.multiply(3);
    EXPECT_EQ(right_hand_side(2, 0).compare(three_times_determinant), 0);
}

TEST(FractionFreeLDLTFactorizationTest, MatchesFlintOnDeterministicSymmetricSamples) {
    std::uint64_t state = 0x8f3f73b5cf1c9adeULL;
    const auto next_value = [&state]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<slong>((state >> 32) % 9) - 4;
    };

    for (size_t dimension = 1; dimension <= 7; ++dimension) {
        for (size_t sample = 0; sample < 40; ++sample) {
            SCOPED_TRACE(::testing::Message() << "dimension=" << dimension << ", sample=" << sample);
            linalg::matrix_int original(dimension, dimension);
            for (size_t row = 0; row < dimension; ++row) {
                for (size_t column = 0; column <= row; ++column) {
                    fmpz_set_si(original(row, column).native_handle(), next_value());
                    original(column, row) = original(row, column);
                }
            }

            linalg::integer expected_determinant;
            fmpz_mat_det(expected_determinant.native_handle(), original.native_handle());
            const bool expected_positive_definite = fmpz_mat_is_spd(original.native_handle()) != 0;
            linalg::matrix_int factored_system(original);
            linalg::fraction_free_ldlt_factorization factorization(dimension);

            ASSERT_EQ(factorization.factorize_inplace(factored_system), expected_determinant.is_zero() ? 0 : 1);
            EXPECT_EQ(factorization.determinant().compare(expected_determinant), 0);
            EXPECT_EQ(factorization.is_positive_definite(), expected_positive_definite);
            if (expected_determinant.is_zero()) continue;

            linalg::matrix_int right_hand_sides(dimension, 2);
            for (size_t row = 0; row < dimension; ++row) {
                for (size_t column = 0; column < 2; ++column) {
                    fmpz_set_si(right_hand_sides(row, column).native_handle(), next_value());
                }
            }
            linalg::matrix_int original_right_hand_sides(right_hand_sides);
            linalg::integer denominator;
            factorization.solve_inplace(right_hand_sides, denominator, factored_system);

            linalg::matrix_int actual(dimension, 2);
            linalg::matrix_int expected(dimension, 2);
            fmpz_mat_mul(actual.native_handle(), original.native_handle(), right_hand_sides.native_handle());
            fmpz_mat_scalar_mul_fmpz(expected.native_handle(), original_right_hand_sides.native_handle(),
                                     denominator.native_handle());
            EXPECT_TRUE(fmpz_mat_equal(actual.native_handle(), expected.native_handle()));
        }
    }
}

TEST(LinearSolverFractionTest, BuildsVectorOnlyForRequestedSuccessfulCandidate) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::zero(); game(0, 1) = fraction::one();
    game(1, 0) = fraction::one();  game(1, 1) = fraction::zero();
    candidate_search::find_candidate_safe finder(game);

    candidate rejected;
    EXPECT_FALSE(finder.find(0b01, 1, rejected, true));
    EXPECT_EQ(rejected.vector.rows(), 0);

    candidate without_vector;
    ASSERT_TRUE(finder.find(0b11, 2, without_vector, false));
    EXPECT_EQ(without_vector.vector.rows(), 0);

    candidate with_vector;
    ASSERT_TRUE(finder.find(0b11, 2, with_vector, true));
    ASSERT_EQ(with_vector.vector.rows(), 2);
    ASSERT_EQ(with_vector.vector.cols(), 1);
    EXPECT_EQ(with_vector.vector(0, 0), fraction(1, 2));
    EXPECT_EQ(with_vector.vector(1, 0), fraction(1, 2));
}
