#include <gtest/gtest.h>
#include <linalg/linear_solver.hpp>

using namespace linalg;

/*
 * Exact linear solver tests for bordered systems used in candidate search.
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

TEST(LinearSolverFractionTest, RejectsNonPositiveSupportVariable) {
    matrix_frc Ab(3, 4);
    Ab(0, 0) = fraction::one();  Ab(0, 1) = fraction::zero(); Ab(0, 2) = fraction::zero(); Ab(0, 3) = fraction::zero();
    Ab(1, 0) = fraction::zero(); Ab(1, 1) = fraction::one();  Ab(1, 2) = fraction::zero(); Ab(1, 3) = fraction::one();
    Ab(2, 0) = fraction::zero(); Ab(2, 1) = fraction::zero(); Ab(2, 2) = fraction::one();  Ab(2, 3) = fraction(-1);

    matrix_frc x;
    EXPECT_FALSE(solve_linear_frc(Ab, x));
}

TEST(LinearSolverFractionTest, RejectsSingularSystem) {
    matrix_frc Ab(2, 3);
    Ab(0, 0) = fraction::one(); Ab(0, 1) = fraction::one(); Ab(0, 2) = fraction::one();
    Ab(1, 0) = fraction::two(); Ab(1, 1) = fraction::two(); Ab(1, 2) = fraction::two();

    matrix_frc x;
    EXPECT_FALSE(solve_linear_frc(Ab, x));
}
