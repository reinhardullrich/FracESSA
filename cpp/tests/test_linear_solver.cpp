#include <gtest/gtest.h>

#include <fracessa/find_candidate_exact.hpp>

/*
 * Exact reduced-Hessian LDL^T tests through the candidate procedure that owns it.
 */

TEST(LinearSolverFractionTest, SimpleCandidate) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::zero(); game(0, 1) = fraction::one();
    game(1, 0) = fraction::one();  game(1, 1) = fraction::zero();

    candidate result;
    candidate_search::find_candidate_exact finder(game);
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
    candidate_search::find_candidate_exact finder(game);
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
    candidate_search::find_candidate_exact finder(game);
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
    candidate_search::find_candidate_exact finder(game);
    EXPECT_FALSE(finder.find(0b11, 2, result, false));
    EXPECT_FALSE(finder.reduced_hessian_is_negative_definite());
}

TEST(LinearSolverFractionTest, RejectsSingularSystem) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::one(); game(0, 1) = fraction::one();
    game(1, 0) = fraction::one(); game(1, 1) = fraction::one();

    candidate result;
    candidate_search::find_candidate_exact finder(game);
    EXPECT_FALSE(finder.find(0b11, 2, result, false));
}
