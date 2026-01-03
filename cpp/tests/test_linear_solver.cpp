#include <gtest/gtest.h>
#include <rational_linalg/linear_solver_fraction.hpp>
#include <rational_linalg/linear_solver_double.hpp>

using namespace rational_linalg;

TEST(LinearSolverFractionTest, SimpleSystem) {
    matrix_fraction Ab(2, 3);
    Ab(0, 0) = fraction::two(); Ab(0, 1) = fraction::one(); Ab(0, 2) = fraction(5);
    Ab(1, 0) = fraction::one(); Ab(1, 1) = fraction::neg_one(); Ab(1, 2) = fraction::one();

    linear_solver_fraction solver(Ab);
    matrix_fraction x;
    EXPECT_TRUE(solver.solve(x));
    EXPECT_EQ(x(0, 0), fraction::two());
    EXPECT_EQ(x(1, 0), fraction::one());
}

TEST(LinearSolverDoubleTest, SimpleSystem) {
    matrix_double Ab(2, 3);
    Ab(0, 0) = 2.0; Ab(0, 1) = 1.0; Ab(0, 2) = 5.0;
    Ab(1, 0) = 1.0; Ab(1, 1) = -1.0; Ab(1, 2) = 1.0;

    linear_solver_double solver(Ab);
    matrix_double x;
    EXPECT_TRUE(solver.solve(x));
    EXPECT_NEAR(x(0, 0), 2.0, 1e-9);
    EXPECT_NEAR(x(1, 0), 1.0, 1e-9);
}
