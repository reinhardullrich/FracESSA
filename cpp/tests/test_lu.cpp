#include <gtest/gtest.h>
#include <linalg/lu_factor_fraction.hpp>

using namespace linalg;

/*
 * Exact LU tests:
 * determinant sign/value, inverse consistency, and solve correctness.
 */

TEST(LUFactorFractionTest, Determinant) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::two(); A(0, 1) = fraction::one();
    A(1, 0) = fraction::one(); A(1, 1) = fraction::two();

    LU_Factorization lu(A);
    EXPECT_EQ(lu.determinant(), fraction(3));
}

TEST(LUFactorFractionTest, Inverse) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::two(); A(0, 1) = fraction::one();
    A(1, 0) = fraction::one(); A(1, 1) = fraction::two();

    LU_Factorization lu(A);
    matrix_frc inv = lu.inverse();
    
    // Check A * A^-1 = I
    matrix_frc I(2, 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            fraction sum = fraction::zero();
            for (size_t k = 0; k < 2; ++k) {
                sum += A(i, k) * inv(k, j);
            }
            if (i == j) EXPECT_EQ(sum, fraction::one());
            else EXPECT_EQ(sum, fraction::zero());
        }
    }
}

TEST(LUFactorFractionTest, Solve) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::two(); A(0, 1) = fraction::one();
    A(1, 0) = fraction::one(); A(1, 1) = fraction::two();
    
    matrix_frc b(2, 1);
    b(0, 0) = fraction(5);
    b(1, 0) = fraction(4);

    LU_Factorization lu(A);
    matrix_frc x = lu.solve(b);
    
    EXPECT_EQ(x(0, 0), fraction::two());
    EXPECT_EQ(x(1, 0), fraction::one());
}
