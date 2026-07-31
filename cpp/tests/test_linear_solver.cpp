#include <gtest/gtest.h>

#include <fracessa/find_candidate_exact.hpp>

/*
 * Exact Gaussian-elimination tests through the candidate procedure that owns it.
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
}

TEST(LinearSolverFractionTest, AllowsNegativePayoffVariable) {
    linalg::matrix_frc game(1, 1);
    game(0, 0) = fraction(-2);

    candidate result;
    candidate_search::find_candidate_exact finder(game);
    ASSERT_TRUE(finder.find(0b1, 1, result, false));
    EXPECT_EQ(result.payoff, fraction(-2));
}

TEST(LinearSolverFractionTest, RejectsNonPositiveSupportVariable) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::one(); game(0, 1) = fraction::zero();
    game(1, 0) = fraction::zero(); game(1, 1) = fraction::zero();

    candidate result;
    candidate_search::find_candidate_exact finder(game);
    EXPECT_FALSE(finder.find(0b11, 2, result, false));
}

TEST(LinearSolverFractionTest, RejectsSingularSystem) {
    linalg::matrix_frc game(2, 2);
    game(0, 0) = fraction::one(); game(0, 1) = fraction::one();
    game(1, 0) = fraction::one(); game(1, 1) = fraction::one();

    candidate result;
    candidate_search::find_candidate_exact finder(game);
    EXPECT_FALSE(finder.find(0b11, 2, result, false));
}
