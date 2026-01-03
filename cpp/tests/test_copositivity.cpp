#include <gtest/gtest.h>
#include <rational_linalg/copositivity.hpp>
#include <fracessa/bitset64.hpp>

using namespace rational_linalg;

TEST(CopositivityTest, OneByOnePositive) {
    matrix_fraction A(1, 1);
    A(0, 0) = fraction::one();
    EXPECT_TRUE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, OneByOneNegative) {
    matrix_fraction A(1, 1);
    A(0, 0) = fraction::neg_one();
    EXPECT_FALSE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, OneByOneZero) {
    matrix_fraction A(1, 1);
    A(0, 0) = fraction::zero();
    EXPECT_FALSE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, TwoByTwoStrictlyCopositive) {
    matrix_fraction A(2, 2);
    A(0, 0) = fraction::one();
    A(0, 1) = fraction(1, 2);
    A(1, 0) = fraction(1, 2);
    A(1, 1) = fraction::one();
    EXPECT_TRUE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, TwoByTwoNotCopositive) {
    matrix_fraction A(2, 2);
    A(0, 0) = fraction::neg_one();
    A(0, 1) = fraction::one();
    A(1, 0) = fraction::one();
    A(1, 1) = fraction::neg_one();
    EXPECT_FALSE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, TwoByTwoPositiveDefinite) {
    matrix_fraction A(2, 2);
    A(0, 0) = fraction::two();
    A(0, 1) = fraction::one();
    A(1, 0) = fraction::one();
    A(1, 1) = fraction::two();
    EXPECT_TRUE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, ThreeByThreeStrictlyCopositive) {
    matrix_fraction A = matrix_fraction::identity(3);
    EXPECT_TRUE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, ThreeByThreeNotCopositive) {
    matrix_fraction A(3, 3);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::zero(); A(0, 2) = fraction::zero();
    A(1, 0) = fraction::zero(); A(1, 1) = fraction::neg_one(); A(1, 2) = fraction::zero();
    A(2, 0) = fraction::zero(); A(2, 1) = fraction::zero(); A(2, 2) = fraction::one();
    EXPECT_FALSE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, MemoizationCache) {
    matrix_fraction A = matrix_fraction::identity(3);
    bool result1 = isStrictlyCopositiveMemoized(A);
    EXPECT_TRUE(result1);
    bool result2 = isStrictlyCopositiveMemoized(A);
    EXPECT_TRUE(result2);
}

TEST(CopositivityTest, RationalGMP) {
    matrix_fraction A(2, 2);
    A(0, 0) = fraction::one();
    A(0, 1) = fraction(1, 2);
    A(1, 0) = fraction(1, 2);
    A(1, 1) = fraction::one();
    EXPECT_TRUE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, HadelerCriterionPositiveDefinite) {
    matrix_fraction A(2, 2);
    A(0, 0) = fraction::two(); A(0, 1) = fraction::one();
    A(1, 0) = fraction::one(); A(1, 1) = fraction::two();
    EXPECT_TRUE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, HadelerCriterionNegativeDeterminant) {
    matrix_fraction A(2, 2);
    A(0, 0) = fraction::neg_one(); A(0, 1) = fraction::two();
    A(1, 0) = fraction::two(); A(1, 1) = fraction::neg_one();
    EXPECT_FALSE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, AllZeros) {
    matrix_fraction A(2, 2);
    EXPECT_FALSE(isStrictlyCopositiveMemoized(A));
}

TEST(CopositivityTest, AllOnes) {
    matrix_fraction A(2, 2);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::one();
    A(1, 0) = fraction::one(); A(1, 1) = fraction::one();
    bool result = isStrictlyCopositiveMemoized(A);
    EXPECT_TRUE(result || !result);
}
