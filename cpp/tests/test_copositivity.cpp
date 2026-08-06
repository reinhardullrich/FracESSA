#include <gtest/gtest.h>
#include <linalg/copositive_fraction.hpp>
#include <fracessa/bitset64.hpp>

using namespace linalg;

/*
 * Copositivity regression tests for exact-rational checker.
 *
 * Cases cover minimal dimensions and representative sign patterns where the
 * Hadeler-based logic should clearly accept or reject.
 */

TEST(CopositivityTest, OneByOnePositive) {
    matrix_frc A(1, 1);
    A(0, 0) = fraction::one();
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, OneByOneNegative) {
    matrix_frc A(1, 1);
    A(0, 0) = fraction::neg_one();
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, OneByOneZero) {
    matrix_frc A(1, 1);
    A(0, 0) = fraction::zero();
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, TwoByTwoStrictlyCopositive) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one();
    A(0, 1) = fraction(1, 2);
    A(1, 0) = fraction(1, 2);
    A(1, 1) = fraction::one();
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, TwoByTwoNotCopositive) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::neg_one();
    A(0, 1) = fraction::one();
    A(1, 0) = fraction::one();
    A(1, 1) = fraction::neg_one();
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, TwoByTwoPositiveDefinite) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::two();
    A(0, 1) = fraction::one();
    A(1, 0) = fraction::one();
    A(1, 1) = fraction::two();
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, ThreeByThreeStrictlyCopositive) {
    matrix_frc A;
    A.set_identity(3);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, ThreeByThreeNotCopositive) {
    matrix_frc A(3, 3);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::zero(); A(0, 2) = fraction::zero();
    A(1, 0) = fraction::zero(); A(1, 1) = fraction::neg_one(); A(1, 2) = fraction::zero();
    A(2, 0) = fraction::zero(); A(2, 1) = fraction::zero(); A(2, 2) = fraction::one();
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, LargeMatrixRejectsBeforeEnumeratingAllSubsets) {
    matrix_frc A(62, 62);
    A(0, 0) = fraction::neg_one();
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, MemoizationCache) {
    matrix_frc A;
    A.set_identity(3);
    bool result1 = is_strictly_copositive(A);
    EXPECT_TRUE(result1);
    bool result2 = is_strictly_copositive(A);
    EXPECT_TRUE(result2);
}

TEST(CopositivityTest, RationalGMP) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one();
    A(0, 1) = fraction(1, 2);
    A(1, 0) = fraction(1, 2);
    A(1, 1) = fraction::one();
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, HadelerCriterionPositiveDefinite) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::two(); A(0, 1) = fraction::one();
    A(1, 0) = fraction::one(); A(1, 1) = fraction::two();
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, HadelerCriterionNegativeDeterminant) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::neg_one(); A(0, 1) = fraction::two();
    A(1, 0) = fraction::two(); A(1, 1) = fraction::neg_one();
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, HadelerCriterionSingularPositiveAdjugate) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one();     A(0, 1) = fraction::neg_one();
    A(1, 0) = fraction::neg_one(); A(1, 1) = fraction::one();
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, AllZeros) {
    matrix_frc A(2, 2);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, AllOnes) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::one();
    A(1, 0) = fraction::one(); A(1, 1) = fraction::one();
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, SharedSignScan) {
    matrix_frc A(4, 4);
    A(0, 0) = fraction::one(); A(0, 1) = fraction::neg_one(); A(0, 2) = fraction::zero(); A(0, 3) = fraction::two();
    A(1, 0) = fraction::neg_one(); A(1, 1) = fraction::two(); A(1, 2) = fraction(-3); A(1, 3) = fraction(4);
    A(2, 0) = fraction::zero(); A(2, 1) = fraction(-3); A(2, 2) = fraction(3); A(2, 3) = fraction(-5);
    A(3, 0) = fraction::two(); A(3, 1) = fraction(4); A(3, 2) = fraction(-5); A(3, 3) = fraction(4);

    const CopositivitySignScan scan = scan_copositivity_signs(A);

    EXPECT_TRUE(scan.all_diagonal_positive);
    EXPECT_EQ(scan.negative_neighbors[0], bs64::single_bit_at_pos(1));
    EXPECT_EQ(scan.negative_neighbors[1], bs64::single_bit_at_pos(0) | bs64::single_bit_at_pos(2));
    EXPECT_EQ(scan.negative_neighbors[2], bs64::single_bit_at_pos(1) | bs64::single_bit_at_pos(3));
    EXPECT_EQ(scan.negative_neighbors[3], bs64::single_bit_at_pos(2));
    EXPECT_EQ(scan.rows_with_negative_off_diagonal, bs64::set_all_n_bits(4));
    EXPECT_EQ(scan.rows_with_positive_off_diagonal,
              bs64::single_bit_at_pos(0) | bs64::single_bit_at_pos(1) | bs64::single_bit_at_pos(3));
    EXPECT_FALSE(scan.all_negative_part_row_sums_positive);
    const linalg::integer expected_sum(4);
    EXPECT_EQ(scan.all_ones_quadratic_value.compare(expected_sum), 0);
}

TEST(CopositivityTest, NegativePartDiagonalDominance) {
    matrix_frc A(4, 4);
    A(0, 0) = fraction(2);  A(0, 1) = fraction(-1); A(0, 2) = fraction(3);  A(0, 3) = fraction::zero();
    A(1, 0) = fraction(-1); A(1, 1) = fraction(4);  A(1, 2) = fraction(-2); A(1, 3) = fraction::one();
    A(2, 0) = fraction(3);  A(2, 1) = fraction(-2); A(2, 2) = fraction(4);  A(2, 3) = fraction(-1);
    A(3, 0) = fraction::zero(); A(3, 1) = fraction::one(); A(3, 2) = fraction(-1); A(3, 3) = fraction(2);

    EXPECT_TRUE(scan_copositivity_signs(A).all_negative_part_row_sums_positive);

    A(0, 0) = fraction::one(); // Row 0 now has the non-strict boundary 1 + (-1) = 0.
    EXPECT_FALSE(scan_copositivity_signs(A).all_negative_part_row_sums_positive);
}

TEST(CopositivityTest, SignScanStopsAtNonPositiveDiagonal) {
    matrix_frc A(2, 2);
    A(0, 0) = fraction::one();
    A(0, 1) = A(1, 0) = fraction::neg_one();
    A(1, 1) = fraction::zero();

    const CopositivitySignScan scan = scan_copositivity_signs(A);

    EXPECT_FALSE(scan.all_diagonal_positive);
    EXPECT_EQ(scan.rows_with_negative_off_diagonal, 0);
    EXPECT_EQ(scan.negative_neighbors[0], 0);
    EXPECT_EQ(scan.negative_neighbors[1], 0);
}

TEST(CopositivityTest, DirectSmallCriteriaCoverBoundaryCases) {
    EXPECT_TRUE(is_strictly_copositive_1x1(fraction::one()));
    EXPECT_FALSE(is_strictly_copositive_1x1(fraction::zero()));

    EXPECT_TRUE(is_strictly_copositive_2x2(fraction::one(), fraction::two(), fraction::one()));
    EXPECT_TRUE(is_strictly_copositive_2x2(fraction::two(), fraction::neg_one(), fraction::two()));
    EXPECT_FALSE(is_strictly_copositive_2x2(fraction::one(), fraction::neg_one(), fraction::one()));
    EXPECT_FALSE(is_strictly_copositive_2x2(fraction::one(), fraction(-2), fraction::one()));

    EXPECT_TRUE(is_strictly_copositive_3x3(
        fraction::one(), fraction::one(), fraction::one(), fraction::one(), fraction::one(), fraction::one()));
    EXPECT_TRUE(is_strictly_copositive_3x3(
        fraction::one(), fraction::two(), fraction::two(), fraction::one(), fraction::two(), fraction::one()));
    EXPECT_FALSE(is_strictly_copositive_3x3(
        fraction::two(), fraction::neg_one(), fraction::neg_one(), fraction::two(), fraction::neg_one(), fraction::two()));
    EXPECT_FALSE(is_strictly_copositive_3x3(
        fraction::one(), fraction(-2), fraction::zero(), fraction::one(), fraction::zero(), fraction::one()));
}
