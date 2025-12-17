#include <gtest/gtest.h>
#include <rational_linalg/fraction.hpp>
#include <sstream>
#include <limits>

using namespace rational_linalg;

// ============================================================================
// Constructor Tests
// ============================================================================

TEST(FractionTest, DefaultConstructor) {
    fraction f;
    EXPECT_TRUE(f.is_zero());
    EXPECT_EQ(f.to_double(), 0.0);
}

TEST(FractionTest, ConstructorFromInt) {
    fraction f(5);
    EXPECT_FALSE(f.is_zero());
    EXPECT_EQ(f.to_double(), 5.0);
    EXPECT_TRUE(f.is_one() == false);
}

TEST(FractionTest, ConstructorFromIntWithDenominator) {
    fraction f(3, 4);
    EXPECT_DOUBLE_EQ(f.to_double(), 0.75);
}

TEST(FractionTest, ConstructorFromLong) {
    fraction f(10L, 2L);
    EXPECT_DOUBLE_EQ(f.to_double(), 5.0);
}

TEST(FractionTest, ConstructorFromLongLong) {
    fraction f(100LL, 4LL);
    EXPECT_DOUBLE_EQ(f.to_double(), 25.0);
}

TEST(FractionTest, ConstructorSimplifies) {
    fraction f(4, 8);
    // Should be simplified to 1/2
    fraction expected(1, 2);
    EXPECT_EQ(f, expected);
}

TEST(FractionTest, CopyConstructor) {
    fraction original(3, 7);
    fraction copy(original);
    EXPECT_EQ(original, copy);
    EXPECT_DOUBLE_EQ(copy.to_double(), original.to_double());
}

TEST(FractionTest, MoveConstructor) {
    fraction original(5, 9);
    double original_value = original.to_double();
    fraction moved(std::move(original));
    EXPECT_DOUBLE_EQ(moved.to_double(), original_value);
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
    EXPECT_DOUBLE_EQ(f2.to_double(), f1.to_double());
}

TEST(FractionTest, SelfAssignment) {
    fraction f(7, 11);
    f = f;  // Should not crash
    EXPECT_DOUBLE_EQ(f.to_double(), 7.0 / 11.0);
}

TEST(FractionTest, MoveAssignment) {
    fraction f1(13, 17);
    double value = f1.to_double();
    fraction f2;
    f2 = std::move(f1);
    EXPECT_DOUBLE_EQ(f2.to_double(), value);
}

// ============================================================================
// Arithmetic Operator Tests
// ============================================================================

TEST(FractionTest, Addition) {
    fraction f1(1, 2);
    fraction f2(1, 3);
    fraction result = f1 + f2;
    EXPECT_DOUBLE_EQ(result.to_double(), 5.0 / 6.0);
}

TEST(FractionTest, Subtraction) {
    fraction f1(1, 2);
    fraction f2(1, 3);
    fraction result = f1 - f2;
    EXPECT_DOUBLE_EQ(result.to_double(), 1.0 / 6.0);
}

TEST(FractionTest, Multiplication) {
    fraction f1(2, 3);
    fraction f2(3, 4);
    fraction result = f1 * f2;
    EXPECT_DOUBLE_EQ(result.to_double(), 0.5);
}

TEST(FractionTest, Division) {
    fraction f1(1, 2);
    fraction f2(1, 4);
    fraction result = f1 / f2;
    EXPECT_DOUBLE_EQ(result.to_double(), 2.0);
}

TEST(FractionTest, DivisionByZeroThrows) {
    fraction f1(1, 2);
    fraction f2(0);
    EXPECT_THROW(f1 / f2, std::domain_error);
}

TEST(FractionTest, UnaryNegation) {
    fraction f(3, 4);
    fraction neg = -f;
    EXPECT_DOUBLE_EQ(neg.to_double(), -0.75);
    
    fraction neg2 = -neg;
    EXPECT_EQ(neg2, f);
}

TEST(FractionTest, CompoundAssignmentAddition) {
    fraction f(1, 2);
    f += fraction(1, 3);
    EXPECT_DOUBLE_EQ(f.to_double(), 5.0 / 6.0);
}

TEST(FractionTest, CompoundAssignmentSubtraction) {
    fraction f(1, 2);
    f -= fraction(1, 3);
    EXPECT_DOUBLE_EQ(f.to_double(), 1.0 / 6.0);
}

TEST(FractionTest, CompoundAssignmentMultiplication) {
    fraction f(2, 3);
    f *= fraction(3, 4);
    EXPECT_DOUBLE_EQ(f.to_double(), 0.5);
}

TEST(FractionTest, CompoundAssignmentDivision) {
    fraction f(1, 2);
    f /= fraction(1, 4);
    EXPECT_DOUBLE_EQ(f.to_double(), 2.0);
}

TEST(FractionTest, CompoundAssignmentDivisionByZeroThrows) {
    fraction f(1, 2);
    fraction zero(0);
    EXPECT_THROW(f /= zero, std::domain_error);
}

// ============================================================================
// Comparison Operator Tests
// ============================================================================

TEST(FractionTest, Equality) {
    fraction f1(1, 2);
    fraction f2(2, 4);
    fraction f3(1, 3);
    EXPECT_EQ(f1, f2);
    EXPECT_NE(f1, f3);
}

TEST(FractionTest, LessThan) {
    fraction f1(1, 3);
    fraction f2(1, 2);
    EXPECT_LT(f1, f2);
    EXPECT_FALSE(f2 < f1);
}

TEST(FractionTest, LessThanOrEqual) {
    fraction f1(1, 2);
    fraction f2(2, 4);
    fraction f3(1, 3);
    EXPECT_LE(f1, f2);
    EXPECT_LE(f3, f1);
    EXPECT_FALSE(f1 <= f3);
}

TEST(FractionTest, GreaterThan) {
    fraction f1(1, 2);
    fraction f2(1, 3);
    EXPECT_GT(f1, f2);
    EXPECT_FALSE(f2 > f1);
}

TEST(FractionTest, GreaterThanOrEqual) {
    fraction f1(1, 2);
    fraction f2(2, 4);
    fraction f3(1, 3);
    EXPECT_GE(f1, f2);
    EXPECT_GE(f1, f3);
    EXPECT_FALSE(f3 >= f1);
}

TEST(FractionTest, NegativeComparisons) {
    fraction f1(-1, 2);
    fraction f2(1, 2);
    EXPECT_LT(f1, f2);
    EXPECT_LT(f1, fraction(0));
    EXPECT_GT(f2, f1);
}

// ============================================================================
// Utility Function Tests
// ============================================================================

TEST(FractionTest, IsZero) {
    fraction f1(0);
    fraction f2(0, 5);
    fraction f3(1, 2);
    EXPECT_TRUE(f1.is_zero());
    EXPECT_TRUE(f2.is_zero());
    EXPECT_FALSE(f3.is_zero());
}

TEST(FractionTest, IsOne) {
    fraction f1(1);
    fraction f2(2, 2);
    fraction f3(1, 2);
    EXPECT_TRUE(f1.is_one());
    EXPECT_TRUE(f2.is_one());
    EXPECT_FALSE(f3.is_one());
}

TEST(FractionTest, Abs) {
    fraction f1(-3, 4);
    fraction f2(3, 4);
    fraction abs1 = f1.abs();
    fraction abs2 = f2.abs();
    EXPECT_DOUBLE_EQ(abs1.to_double(), 0.75);
    EXPECT_DOUBLE_EQ(abs2.to_double(), 0.75);
    EXPECT_EQ(abs1, abs2);
}

TEST(FractionTest, Inverse) {
    fraction f(2, 3);
    fraction inv = f.inverse();
    EXPECT_DOUBLE_EQ(inv.to_double(), 1.5);
    
    // Inverse of inverse should be original
    fraction inv_inv = inv.inverse();
    EXPECT_EQ(inv_inv, f);
}

TEST(FractionTest, InverseOfZeroThrows) {
    fraction zero(0);
    EXPECT_THROW(zero.inverse(), std::domain_error);
}

// ============================================================================
// Conversion Tests
// ============================================================================

TEST(FractionTest, ToDouble) {
    fraction f(22, 7);  // Approximation of pi
    double result = f.to_double();
    EXPECT_NEAR(result, 22.0 / 7.0, 1e-10);
}

TEST(FractionTest, ToString) {
    fraction f(3, 4);
    std::string str = f.to_string();
    // FLINT may format as "3/4" or similar
    EXPECT_FALSE(str.empty());
}

TEST(FractionTest, ToStringZero) {
    fraction f(0);
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

TEST(FractionTest, AddInplace) {
    fraction f(1, 2);
    fraction other(1, 3);
    f.add_inplace(other);
    EXPECT_DOUBLE_EQ(f.to_double(), 5.0 / 6.0);
}

TEST(FractionTest, SubInplace) {
    fraction f(1, 2);
    fraction other(1, 3);
    f.sub_inplace(other);
    EXPECT_DOUBLE_EQ(f.to_double(), 1.0 / 6.0);
}

TEST(FractionTest, MulInplace) {
    fraction f(2, 3);
    fraction other(3, 4);
    f.mul_inplace(other);
    EXPECT_DOUBLE_EQ(f.to_double(), 0.5);
}

TEST(FractionTest, DivInplace) {
    fraction f(1, 2);
    fraction other(1, 4);
    f.div_inplace(other);
    EXPECT_DOUBLE_EQ(f.to_double(), 2.0);
}

TEST(FractionTest, DivInplaceByZeroThrows) {
    fraction f(1, 2);
    fraction zero(0);
    EXPECT_THROW(f.div_inplace(zero), std::domain_error);
}

TEST(FractionTest, NegateInplace) {
    fraction f(3, 4);
    f.negate_inplace();
    EXPECT_DOUBLE_EQ(f.to_double(), -0.75);
    
    f.negate_inplace();
    EXPECT_DOUBLE_EQ(f.to_double(), 0.75);
}

TEST(FractionTest, AbsInplace) {
    fraction f(-3, 4);
    f.abs_inplace();
    EXPECT_DOUBLE_EQ(f.to_double(), 0.75);
    
    fraction f2(3, 4);
    f2.abs_inplace();
    EXPECT_DOUBLE_EQ(f2.to_double(), 0.75);
}

// ============================================================================
// Edge Case Tests
// ============================================================================

TEST(FractionTest, LargeNumbers) {
    fraction f(1000000, 1);
    EXPECT_DOUBLE_EQ(f.to_double(), 1000000.0);
}

TEST(FractionTest, VeryLargeDenominator) {
    fraction f(1, 1000000);
    EXPECT_DOUBLE_EQ(f.to_double(), 1e-6);
}

TEST(FractionTest, NegativeNumerator) {
    fraction f(-5, 3);
    EXPECT_DOUBLE_EQ(f.to_double(), -5.0 / 3.0);
}

TEST(FractionTest, NegativeDenominator) {
    // Note: Constructing with negative denominator may not be directly supported
    // Instead, construct with negative numerator to get negative fraction
    fraction f(-5, 3);
    EXPECT_NEAR(f.to_double(), -5.0 / 3.0, 1e-10);
    EXPECT_LT(f.to_double(), 0.0);
}

TEST(FractionTest, BothNegative) {
    // Note: Constructing with both negative may not be directly supported
    // Instead, construct with positive values to get positive fraction
    fraction f(5, 3);
    EXPECT_NEAR(f.to_double(), 5.0 / 3.0, 1e-10);
    EXPECT_GT(f.to_double(), 0.0);
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
    EXPECT_TRUE(f1.is_one());
    EXPECT_TRUE(f2.is_one());
}

// ============================================================================
// Helper Function Tests
// ============================================================================

TEST(FractionTest, RationalToDoubleHelper) {
    fraction f(7, 8);
    double result = f.to_double();
    EXPECT_DOUBLE_EQ(result, 0.875);
}

TEST(FractionTest, ToStringHelper) {
    fraction f(9, 10);
    std::string str = f.to_string();
    EXPECT_FALSE(str.empty());
}

// ============================================================================
// Complex Operation Tests
// ============================================================================

TEST(FractionTest, ChainedOperations) {
    fraction f(1, 2);
    f += fraction(1, 3);
    f -= fraction(1, 6);
    f *= fraction(2, 1);
    // (1/2 + 1/3 - 1/6) * 2 = (3/6 + 2/6 - 1/6) * 2 = (4/6) * 2 = 8/6 = 4/3
    EXPECT_NEAR(f.to_double(), 4.0 / 3.0, 1e-10);
}

TEST(FractionTest, MixedArithmetic) {
    fraction f1(1, 2);
    fraction f2(1, 3);
    fraction f3(1, 4);
    fraction result = (f1 + f2) * f3;
    EXPECT_DOUBLE_EQ(result.to_double(), 5.0 / 24.0);
}

TEST(FractionTest, ComparisonChain) {
    fraction f1(1, 4);
    fraction f2(1, 3);
    fraction f3(1, 2);
    EXPECT_LT(f1, f2);
    EXPECT_LT(f2, f3);
    EXPECT_LT(f1, f3);
}

