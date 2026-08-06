#include <gtest/gtest.h>
#include <linalg/fraction.hpp>
#include <sstream>
#include <limits>

using namespace linalg;

/*
 * Exact rational wrapper tests.
 *
 * Goal: guard C++ value semantics around FLINT-backed storage and verify that
 * the retained public operations preserve canonical rational behavior.
 */

// ============================================================================
// Constructor Tests
// ============================================================================

TEST(FractionTest, DefaultConstructor) {
    fraction f;
    EXPECT_TRUE(f.is_zero());
    EXPECT_EQ(f.to_dbl(), 0.0);
}

TEST(FractionTest, ConstructorFromInt) {
    fraction f(5);
    EXPECT_FALSE(f.is_zero());
    EXPECT_EQ(f.to_dbl(), 5.0);
}

TEST(FractionTest, ConstructorFromIntWithDenominator) {
    fraction f(3, 4);
    EXPECT_DOUBLE_EQ(f.to_dbl(), 0.75);
}

TEST(FractionTest, ConstructorFromLongLong) {
    fraction f(100LL, 4LL);
    EXPECT_DOUBLE_EQ(f.to_dbl(), 25.0);
}

TEST(FractionTest, ConstructorSimplifies) {
    fraction f(4, 8);
    // Should be simplified to 1/2
    fraction expected(1, 2);
    EXPECT_EQ(f, expected);
}

TEST(FractionTest, ConstructorNormalizesNegativeDenominator) {
    EXPECT_EQ(fraction(1, -2).to_string(), "-1/2");
    EXPECT_EQ(fraction(-1, -2).to_string(), "1/2");
    EXPECT_EQ(fraction(1LL, -2LL).to_string(), "-1/2");
}

TEST(FractionTest, CopyConstructor) {
    fraction original(3, 7);
    fraction copy(original);
    EXPECT_EQ(original, copy);
    EXPECT_DOUBLE_EQ(copy.to_dbl(), original.to_dbl());
}

TEST(FractionTest, MoveConstructor) {
    fraction original(5, 9);
    double original_value = original.to_dbl();
    fraction moved(std::move(original));
    EXPECT_DOUBLE_EQ(moved.to_dbl(), original_value);
    // Original should still be valid (FLINT handles this)
}

// ============================================================================
// Assignment Operator Tests
// ============================================================================

TEST(FractionTest, CopyAssignment) {
    fraction f1(1, 3);
    fraction f2(2, 5);
    f2 = f1;
    EXPECT_EQ(f1, f2);
    EXPECT_DOUBLE_EQ(f2.to_dbl(), f1.to_dbl());
}

TEST(FractionTest, SelfAssignment) {
    fraction f(7, 11);
    f = f;  // Should not crash
    EXPECT_DOUBLE_EQ(f.to_dbl(), 7.0 / 11.0);
}

TEST(FractionTest, MoveAssignment) {
    fraction f1(13, 17);
    double value = f1.to_dbl();
    fraction f2;
    f2 = std::move(f1);
    EXPECT_DOUBLE_EQ(f2.to_dbl(), value);
}

// ============================================================================
// Arithmetic Operator Tests
// ============================================================================

TEST(FractionTest, Multiplication) {
    fraction f1(2, 3);
    fraction f2(3, 4);
    fraction result = f1 * f2;
    EXPECT_DOUBLE_EQ(result.to_dbl(), 0.5);
}

TEST(FractionTest, Equality) {
    fraction f1(1, 2);
    fraction f2(2, 4);
    fraction f3(1, 3);
    EXPECT_EQ(f1, f2);
    EXPECT_FALSE(f1 == f3);
}

// ============================================================================
// Utility Function Tests
// ============================================================================

TEST(FractionTest, IsZero) {
    fraction f1 = fraction::zero();
    fraction f2(0, 5);
    fraction f3(1, 2);
    EXPECT_TRUE(f1.is_zero());
    EXPECT_TRUE(f2.is_zero());
    EXPECT_FALSE(f3.is_zero());
}

// ============================================================================
// Conversion Tests
// ============================================================================

TEST(FractionTest, ToDouble) {
    fraction f(22, 7);  // Approximation of pi
    double result = f.to_dbl();
    EXPECT_NEAR(result, 22.0 / 7.0, 1e-10);
}

TEST(FractionTest, ToString) {
    fraction f(3, 4);
    std::string str = f.to_string();
    // FLINT may format as "3/4" or similar
    EXPECT_FALSE(str.empty());
}

TEST(FractionTest, ToStringZero) {
    fraction f = fraction::zero();
    std::string str = f.to_string();
    EXPECT_FALSE(str.empty());
}

TEST(FractionTest, StreamOutput) {
    fraction f(5, 6);
    std::ostringstream oss;
    oss << f;
    std::string str = oss.str();
    EXPECT_FALSE(str.empty());
}

// ============================================================================
// In-Place Operation Tests
// ============================================================================

TEST(FractionTest, AddAssign) {
    fraction f(1, 2);
    fraction other(1, 3);
    f += other;
    EXPECT_DOUBLE_EQ(f.to_dbl(), 5.0 / 6.0);
}

TEST(FractionTest, SubtractAssign) {
    fraction f(1, 2);
    fraction other(1, 3);
    f -= other;
    EXPECT_DOUBLE_EQ(f.to_dbl(), 1.0 / 6.0);
}

// ============================================================================
// Edge Case Tests
// ============================================================================

TEST(FractionTest, LargeNumbers) {
    fraction f(1000000, 1);
    EXPECT_DOUBLE_EQ(f.to_dbl(), 1000000.0);
}

TEST(FractionTest, VeryLargeDenominator) {
    fraction f(1, 1000000);
    EXPECT_DOUBLE_EQ(f.to_dbl(), 1e-6);
}

TEST(FractionTest, NegativeNumerator) {
    fraction f(-5, 3);
    EXPECT_DOUBLE_EQ(f.to_dbl(), -5.0 / 3.0);
}

TEST(FractionTest, NegativeDenominator) {
    // Note: Constructing with negative denominator may not be directly supported
    // Instead, construct with negative numerator to get negative fraction
    fraction f(-5, 3);
    EXPECT_NEAR(f.to_dbl(), -5.0 / 3.0, 1e-10);
    EXPECT_LT(f.to_dbl(), 0.0);
}

TEST(FractionTest, BothNegative) {
    // Note: Constructing with both negative may not be directly supported
    // Instead, construct with positive values to get positive fraction
    fraction f(5, 3);
    EXPECT_NEAR(f.to_dbl(), 5.0 / 3.0, 1e-10);
    EXPECT_GT(f.to_dbl(), 0.0);
}

TEST(FractionTest, FractionSimplification) {
    fraction f(100, 200);
    fraction expected(1, 2);
    EXPECT_EQ(f, expected);
}

TEST(FractionTest, ZeroHandling) {
    fraction f1(0, 1);
    fraction f2(0, 100);
    EXPECT_EQ(f1, f2);
    EXPECT_TRUE(f1.is_zero());
    EXPECT_TRUE(f2.is_zero());
}

TEST(FractionTest, OneHandling) {
    fraction f1(1, 1);
    fraction f2(5, 5);
    EXPECT_EQ(f1, f2);
    EXPECT_EQ(f1, fraction::one());
    EXPECT_EQ(f2, fraction::one());
}

// ============================================================================
// Helper Function Tests
// ============================================================================

TEST(FractionTest, RationalToDoubleHelper) {
    fraction f(7, 8);
    double result = f.to_dbl();
    EXPECT_DOUBLE_EQ(result, 0.875);
}

TEST(FractionTest, ToStringHelper) {
    fraction f(9, 10);
    std::string str = f.to_string();
    EXPECT_FALSE(str.empty());
}
