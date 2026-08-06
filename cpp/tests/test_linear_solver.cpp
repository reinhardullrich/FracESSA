#include <gtest/gtest.h>

#include <fracessa/find_candidate_safe.hpp>

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
