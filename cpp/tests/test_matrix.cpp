#include <gtest/gtest.h>
#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_double.hpp>

using namespace linalg;

/*
 * Matrix container and factory tests.
 *
 * These keep basic storage/indexing/factory invariants stable and exercise
 * positive-definiteness routines on small, hand-checkable examples.
 */

TEST(MatrixFractionTest, BasicOperations) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::two();
    A(1, 0) = fraction(3); A(1, 1) = fraction(4);

    EXPECT_EQ(A.rows(), 2);
    EXPECT_EQ(A.cols(), 2);
    EXPECT_EQ(A(0, 0), fraction::one());
    EXPECT_EQ(A(1, 1), fraction(4));
}

TEST(MatrixDoubleTest, BasicOperations) {
    matrix_dbl A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 2.0;
    A(1, 0) = 3.0; A(1, 1) = 4.0;

    EXPECT_EQ(A.rows(), 2);
    EXPECT_EQ(A.cols(), 2);
    EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A(1, 1), 4.0);
}

TEST(MatrixFractionTest, FactoryFunctions) {
    std::vector<fraction> vals = {fraction::one(), fraction::two(), fraction(3)};
    matrix_frc S = create_symmetric(2, vals);
    EXPECT_EQ(S(0, 0), fraction::one());
    EXPECT_EQ(S(0, 1), fraction::two());
    EXPECT_EQ(S(1, 0), fraction::two());
    EXPECT_EQ(S(1, 1), fraction(3));
}

TEST(MatrixPositiveDefiniteTest, Fraction) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::two(); A(0, 1) = fraction::neg_one();
    A(1, 0) = fraction::neg_one(); A(1, 1) = fraction::two();
    EXPECT_TRUE(A.is_positive_definite());

    matrix_frc B(2, 2);
    B(0, 0) = fraction::one(); B(0, 1) = fraction::two();
    B(1, 0) = fraction::two(); B(1, 1) = fraction::one();
    EXPECT_FALSE(B.is_positive_definite());
}

TEST(MatrixPositiveDefiniteTest, Double) {
    matrix_dbl A(2, 2);
    A(0, 0) = 2.0; A(0, 1) = -1.0;
    A(1, 0) = -1.0; A(1, 1) = 2.0;
    EXPECT_TRUE(A.is_positive_definite());
}

TEST(MatrixPositiveDefiniteTest, DoubleRejectsIndefinite) {
    matrix_dbl B(2, 2);
    B(0, 0) = 1.0; B(0, 1) = 2.0;
    B(1, 0) = 2.0; B(1, 1) = 1.0;
    EXPECT_FALSE(B.is_positive_definite());
}

TEST(MatrixPositiveDefiniteTest, DoubleRejectsNonSquare) {
    matrix_dbl C(2, 3);
    C(0, 0) = 1.0; C(0, 1) = 0.0; C(0, 2) = 0.0;
    C(1, 0) = 0.0; C(1, 1) = 1.0; C(1, 2) = 0.0;
    EXPECT_FALSE(C.is_positive_definite());
}
