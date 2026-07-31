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
                sum.addmul(A(i, k), inv(k, j));
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

TEST(LUFactorFractionTest, PivotAfterEarlierElimination) {
    matrix_frc A(3, 3);
    A(0, 0) = fraction::one();  A(0, 1) = fraction::one(); A(0, 2) = fraction::zero();
    A(1, 0) = fraction::one();  A(1, 1) = fraction::one(); A(1, 2) = fraction::one();
    A(2, 0) = fraction::zero(); A(2, 1) = fraction::one(); A(2, 2) = fraction::one();

    matrix_frc b(3, 1);
    b(0, 0) = fraction(3);
    b(1, 0) = fraction(6);
    b(2, 0) = fraction(5);

    LU_Factorization lu(A);
    EXPECT_FALSE(lu.isSingular());
    EXPECT_EQ(lu.determinant(), fraction::neg_one());

    const matrix_frc x = lu.solve(b);
    EXPECT_EQ(x(0, 0), fraction::one());
    EXPECT_EQ(x(1, 0), fraction::two());
    EXPECT_EQ(x(2, 0), fraction(3));

    const matrix_frc inverse = lu.inverse();
    for (size_t row = 0; row < 3; ++row) {
        for (size_t column = 0; column < 3; ++column) {
            fraction product = fraction::zero();
            for (size_t k = 0; k < 3; ++k) {
                product.addmul(A(row, k), inverse(k, column));
            }
            EXPECT_EQ(product, row == column ? fraction::one() : fraction::zero());
        }
    }
}

TEST(LUFactorFractionTest, SingularDeterminantAndFlag) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::two();
    A(1, 0) = fraction::two(); A(1, 1) = fraction(4); // row 2 = 2 * row 1

    LU_Factorization lu(A);
    EXPECT_TRUE(lu.isSingular());
    EXPECT_EQ(lu.determinant(), fraction::zero());
}

TEST(LUFactorFractionTest, SingularInverseThrows) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::two();
    A(1, 0) = fraction::two(); A(1, 1) = fraction(4);

    LU_Factorization lu(A);
    EXPECT_THROW(lu.inverse(), std::runtime_error);
}

TEST(LUFactorFractionTest, SingularSolveThrows) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::two();
    A(1, 0) = fraction::two(); A(1, 1) = fraction(4);

    matrix_frc b(2, 1);
    b(0, 0) = fraction::one();
    b(1, 0) = fraction::two();

    LU_Factorization lu(A);
    EXPECT_THROW(lu.solve(b), std::domain_error);
}
