#include <gtest/gtest.h>
#include <fracessa/find_candidate_unsafe.hpp>
#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_double.hpp>

using namespace linalg;
using namespace candidate_search;

/*
 * Matrix container and factory tests.
 *
 * These keep basic storage/indexing/factory invariants stable and exercise
 * exact positive-definiteness on small, hand-checkable examples.
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

TEST(FindCandidateUnsafeTest, SetsInputWarningForEveryRiskyConversion) {
    const auto inspect = [](const fraction& diagonal_0, const fraction& off_diagonal, const fraction& diagonal_1) {
        matrix_frc A(2, 2);
        A(0, 0) = diagonal_0;   A(0, 1) = off_diagonal;
        A(1, 0) = off_diagonal; A(1, 1) = diagonal_1;
        find_candidate_unsafe unsafe(A);
        unsafe.convert_game_matrix();
        return unsafe.input_warnings();
    };

    EXPECT_FALSE(inspect(fraction::zero(), fraction::one(), fraction::zero()));

    EXPECT_TRUE(inspect(fraction::zero(), fraction("1/100000000000000000000"), fraction::zero()));

    EXPECT_TRUE(inspect(fraction("100000000000000000000"), fraction("100000000000000000001"),
                        fraction("100000000000000000000")));

    fraction underflow = fraction::one();
    for (size_t i = 0; i < 1075; ++i) underflow.div_inplace(fraction::two());
    EXPECT_TRUE(inspect(fraction::zero(), underflow, fraction::zero()));

    fraction overflow = fraction::one();
    for (size_t i = 0; i < 1024; ++i) overflow.mul_inplace(fraction::two());
    EXPECT_TRUE(inspect(fraction::zero(), overflow, fraction::zero()));

    fraction subnormal = fraction::one();
    for (size_t i = 0; i < 1023; ++i) subnormal.div_inplace(fraction::two());
    EXPECT_TRUE(inspect(fraction::zero(), subnormal, fraction::zero()));

    fraction large = fraction::one();
    for (size_t i = 0; i < 60; ++i) large.mul_inplace(fraction::two());
    EXPECT_TRUE(inspect(fraction::one(), large, fraction::one()));
}
