#include <gtest/gtest.h>
#include <linalg/linear_solver.hpp>

using namespace linalg;

/*
 * Linear solver unit tests for the bordered systems used in candidate search.
 * Both exact and double paths should recover the same known solution.
 */

TEST(LinearSolverFractionTest, SimpleSystem) {
    matrix_frc Ab(2, 3);
    Ab(0, 0) = fraction::two(); Ab(0, 1) = fraction::one(); Ab(0, 2) = fraction(5);
    Ab(1, 0) = fraction::one(); Ab(1, 1) = fraction::neg_one(); Ab(1, 2) = fraction::one();

    matrix_frc x;
    EXPECT_TRUE(solve_linear_frc(Ab, x));
    EXPECT_EQ(x(0, 0), fraction::two());
    EXPECT_EQ(x(1, 0), fraction::one());
}

TEST(LinearSolverDoubleTest, SimpleSystem) {
    matrix_dbl Ab(2, 3);
    Ab(0, 0) = 2.0; Ab(0, 1) = 1.0; Ab(0, 2) = 5.0;
    Ab(1, 0) = 1.0; Ab(1, 1) = -1.0; Ab(1, 2) = 1.0;

    matrix_dbl x;
    EXPECT_TRUE(solve_linear_dbl(Ab, x));
    EXPECT_NEAR(x(0, 0), 2.0, 1e-9);
    EXPECT_NEAR(x(1, 0), 1.0, 1e-9);
}

TEST(LinearSolverDoubleTest, AllowsNegativePayoffVariable) {
    // x0, x1 are support probabilities (>0), x2 is payoff (may be negative).
    matrix_dbl Ab(3, 4);
    Ab(0, 0) = 1.0; Ab(0, 1) = 0.0; Ab(0, 2) = 0.0; Ab(0, 3) = 0.25;
    Ab(1, 0) = 0.0; Ab(1, 1) = 1.0; Ab(1, 2) = 0.0; Ab(1, 3) = 0.75;
    Ab(2, 0) = 0.0; Ab(2, 1) = 0.0; Ab(2, 2) = 1.0; Ab(2, 3) = -2.0;

    matrix_dbl x;
    ASSERT_TRUE(solve_linear_dbl(Ab, x));
    EXPECT_NEAR(x(0, 0), 0.25, 1e-12);
    EXPECT_NEAR(x(1, 0), 0.75, 1e-12);
    EXPECT_NEAR(x(2, 0), -2.0, 1e-12);
}

TEST(LinearSolverFractionTest, AllowsNegativePayoffVariable) {
    matrix_frc Ab(3, 4);
    Ab(0, 0) = fraction::one();  Ab(0, 1) = fraction::zero(); Ab(0, 2) = fraction::zero(); Ab(0, 3) = fraction(1, 4);
    Ab(1, 0) = fraction::zero(); Ab(1, 1) = fraction::one();  Ab(1, 2) = fraction::zero(); Ab(1, 3) = fraction(3, 4);
    Ab(2, 0) = fraction::zero(); Ab(2, 1) = fraction::zero(); Ab(2, 2) = fraction::one();  Ab(2, 3) = fraction(-2);

    matrix_frc x;
    ASSERT_TRUE(solve_linear_frc(Ab, x));
    EXPECT_EQ(x(0, 0), fraction(1, 4));
    EXPECT_EQ(x(1, 0), fraction(3, 4));
    EXPECT_EQ(x(2, 0), fraction(-2));
}

TEST(LinearSolverDoubleTest, RejectsNonPositiveSupportVariable) {
    matrix_dbl Ab(3, 4);
    Ab(0, 0) = 1.0; Ab(0, 1) = 0.0; Ab(0, 2) = 0.0; Ab(0, 3) = -0.1;
    Ab(1, 0) = 0.0; Ab(1, 1) = 1.0; Ab(1, 2) = 0.0; Ab(1, 3) = 1.1;
    Ab(2, 0) = 0.0; Ab(2, 1) = 0.0; Ab(2, 2) = 1.0; Ab(2, 3) = 0.5;

    matrix_dbl x;
    EXPECT_FALSE(solve_linear_dbl(Ab, x));
}

TEST(LinearSolverFractionTest, RejectsNonPositiveSupportVariable) {
    matrix_frc Ab(3, 4);
    Ab(0, 0) = fraction::one();  Ab(0, 1) = fraction::zero(); Ab(0, 2) = fraction::zero(); Ab(0, 3) = fraction::zero();
    Ab(1, 0) = fraction::zero(); Ab(1, 1) = fraction::one();  Ab(1, 2) = fraction::zero(); Ab(1, 3) = fraction::one();
    Ab(2, 0) = fraction::zero(); Ab(2, 1) = fraction::zero(); Ab(2, 2) = fraction::one();  Ab(2, 3) = fraction(-1);

    matrix_frc x;
    EXPECT_FALSE(solve_linear_frc(Ab, x));
}

TEST(LinearSolverDoubleTest, RejectsSingularSystem) {
    matrix_dbl Ab(2, 3);
    Ab(0, 0) = 1.0; Ab(0, 1) = 1.0; Ab(0, 2) = 1.0;
    Ab(1, 0) = 2.0; Ab(1, 1) = 2.0; Ab(1, 2) = 2.0;

    matrix_dbl x;
    EXPECT_FALSE(solve_linear_dbl(Ab, x));
}

TEST(LinearSolverFractionTest, RejectsSingularSystem) {
    matrix_frc Ab(2, 3);
    Ab(0, 0) = fraction::one(); Ab(0, 1) = fraction::one(); Ab(0, 2) = fraction::one();
    Ab(1, 0) = fraction::two(); Ab(1, 1) = fraction::two(); Ab(1, 2) = fraction::two();

    matrix_frc x;
    EXPECT_FALSE(solve_linear_frc(Ab, x));
}
