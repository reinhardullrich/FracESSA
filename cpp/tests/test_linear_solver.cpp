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
