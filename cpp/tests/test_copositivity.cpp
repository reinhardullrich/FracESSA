#include <gtest/gtest.h>
#include <linalg/copositive_integer.hpp>

using namespace linalg;

/*
 * Copositivity regression tests for exact-integer checker.
 *
 * Cases cover minimal dimensions and representative sign patterns where the
 * Hadeler-based logic should clearly accept or reject.
 */

TEST(CopositivityTest, OneByOnePositive) {
    matrix_int A(1, 1);
    A(0, 0) = integer(1);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, OneByOneNegative) {
    matrix_int A(1, 1);
    A(0, 0) = integer(-1);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, OneByOneZero) {
    matrix_int A(1, 1);
    A(0, 0) = integer(0);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, TwoByTwoStrictlyCopositive) {
    matrix_int A(2, 2);
    A(0, 0) = integer(2);
    A(0, 1) = integer(1);
    A(1, 0) = integer(1);
    A(1, 1) = integer(2);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, TwoByTwoNotCopositive) {
    matrix_int A(2, 2);
    A(0, 0) = integer(-1);
    A(0, 1) = integer(1);
    A(1, 0) = integer(1);
    A(1, 1) = integer(-1);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, TwoByTwoPositiveDefinite) {
    matrix_int A(2, 2);
    A(0, 0) = integer(2);
    A(0, 1) = integer(1);
    A(1, 0) = integer(1);
    A(1, 1) = integer(2);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, ThreeByThreeStrictlyCopositive) {
    matrix_int A;
    A.set_identity(3);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, ThreeByThreeNotCopositive) {
    matrix_int A(3, 3);
    A(0, 0) = integer(1); A(0, 1) = integer(0); A(0, 2) = integer(0);
    A(1, 0) = integer(0); A(1, 1) = integer(-1); A(1, 2) = integer(0);
    A(2, 0) = integer(0); A(2, 1) = integer(0); A(2, 2) = integer(1);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, LargeMatrixRejectsImmediately) {
    matrix_int A(62, 62);
    A(0, 0) = integer(-1);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, RepeatedCallsAreDeterministic) {
    matrix_int A;
    A.set_identity(3);
    bool result1 = is_strictly_copositive(A);
    EXPECT_TRUE(result1);
    bool result2 = is_strictly_copositive(A);
    EXPECT_TRUE(result2);
}

TEST(CopositivityTest, SmallCriterionPositiveDefinite) {
    matrix_int A(2, 2);
    A(0, 0) = integer(2); A(0, 1) = integer(1);
    A(1, 0) = integer(1); A(1, 1) = integer(2);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, SmallCriterionNegativeDeterminant) {
    matrix_int A(2, 2);
    A(0, 0) = integer(-1); A(0, 1) = integer(2);
    A(1, 0) = integer(2); A(1, 1) = integer(-1);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, SmallCriterionSingularPositiveAdjugate) {
    matrix_int A(2, 2);
    A(0, 0) = integer(1);     A(0, 1) = integer(-1);
    A(1, 0) = integer(-1); A(1, 1) = integer(1);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, FourByFourIdentityPasses) {
    matrix_int A;
    A.set_identity(4);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, FourByFourNegativeDeterminantUsesOneRightHandSide) {
    matrix_int A(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) A(row, column) = integer(row == column ? 5 : -2);
    }
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, FourByFourPositiveOffDiagonalPasses) {
    matrix_int A(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) A(row, column) = integer(row == column ? 1 : 2);
    }
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, FourByFourSingularUsesPositiveNullspace) {
    matrix_int A(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) A(row, column) = integer(row == column ? 3 : -1);
    }
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, FourByFourMixedSignKernelPasses) {
    constexpr slong v[4] = {1, -1, 1, -1};
    matrix_int A(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) {
            A(row, column) = integer((row == column ? 4 : 0) - v[row] * v[column]);
        }
    }
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, FourByFourAllOnesPasses) {
    matrix_int A(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) A(row, column) = integer(1);
    }
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, ArbitraryPrecisionIntegerHadelerBranches) {
    integer big;
    ASSERT_EQ(fmpz_set_str(big.native_handle(), "123456789012345678901234567890", 10), 0);

    matrix_int negative_determinant(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) {
            negative_determinant(row, column) = big;
            negative_determinant(row, column).multiply(row == column ? 5 : 2);
            if (row != column) negative_determinant(row, column).negate();
        }
    }
    EXPECT_FALSE(is_strictly_copositive(negative_determinant));

    matrix_int singular(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) {
            singular(row, column) = big;
            singular(row, column).multiply(row == column ? 3 : 1);
            if (row != column) singular(row, column).negate();
        }
    }
    EXPECT_FALSE(is_strictly_copositive(singular));
}

TEST(CopositivityTest, AllZeros) {
    matrix_int A(2, 2);
    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, AllOnes) {
    matrix_int A(2, 2);
    A(0, 0) = integer(1); A(0, 1) = integer(1);
    A(1, 0) = integer(1); A(1, 1) = integer(1);
    EXPECT_TRUE(is_strictly_copositive(A));
}

TEST(CopositivityTest, SharedSignScan) {
    matrix_int A(4, 4);
    A(0, 0) = integer(1); A(0, 1) = integer(-1); A(0, 2) = integer(0); A(0, 3) = integer(2);
    A(1, 0) = integer(-1); A(1, 1) = integer(2); A(1, 2) = integer(-3); A(1, 3) = integer(4);
    A(2, 0) = integer(0); A(2, 1) = integer(-3); A(2, 2) = integer(3); A(2, 3) = integer(-5);
    A(3, 0) = integer(2); A(3, 1) = integer(4); A(3, 2) = integer(-5); A(3, 3) = integer(4);

    const CopositivitySignScan scan = scan_copositivity_signs(A);

    EXPECT_TRUE(scan.all_diagonal_positive);
    EXPECT_TRUE(scan.has_negative_off_diagonal);
    EXPECT_EQ(scan.negative_neighbors[0], bs64::single_bit_at_pos(1));
    EXPECT_EQ(scan.negative_neighbors[1], bs64::single_bit_at_pos(0) | bs64::single_bit_at_pos(2));
    EXPECT_EQ(scan.negative_neighbors[2], bs64::single_bit_at_pos(1) | bs64::single_bit_at_pos(3));
    EXPECT_EQ(scan.negative_neighbors[3], bs64::single_bit_at_pos(2));
    EXPECT_FALSE(scan.all_negative_part_row_sums_positive);
    const linalg::integer expected_sum(4);
    EXPECT_EQ(scan.all_ones_quadratic_value.compare(expected_sum), 0);
}

TEST(CopositivityTest, NegativePartDiagonalDominance) {
    matrix_int A(4, 4);
    A(0, 0) = integer(2);  A(0, 1) = integer(-1); A(0, 2) = integer(3);  A(0, 3) = integer(0);
    A(1, 0) = integer(-1); A(1, 1) = integer(4);  A(1, 2) = integer(-2); A(1, 3) = integer(1);
    A(2, 0) = integer(3);  A(2, 1) = integer(-2); A(2, 2) = integer(4);  A(2, 3) = integer(-1);
    A(3, 0) = integer(0); A(3, 1) = integer(1); A(3, 2) = integer(-1); A(3, 3) = integer(2);

    EXPECT_TRUE(scan_copositivity_signs(A).all_negative_part_row_sums_positive);

    A(0, 0) = integer(1); // Row 0 now has the non-strict boundary 1 + (-1) = 0.
    EXPECT_FALSE(scan_copositivity_signs(A).all_negative_part_row_sums_positive);
}

TEST(CopositivityTest, SignScanStopsAtNonPositiveDiagonal) {
    matrix_int A(2, 2);
    A(0, 0) = integer(1);
    A(0, 1) = A(1, 0) = integer(-1);
    A(1, 1) = integer(0);

    const CopositivitySignScan scan = scan_copositivity_signs(A);

    EXPECT_FALSE(scan.all_diagonal_positive);
    EXPECT_FALSE(scan.has_negative_off_diagonal);
    EXPECT_EQ(scan.negative_neighbors[0], 0U);
    EXPECT_EQ(scan.negative_neighbors[1], 0U);
}

TEST(CopositivityTest, ConnectedNegativeGraphChecksWholeMatrix) {
    matrix_int A(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) A(row, column) = integer(row == column ? 5 : -1);
    }

    const CopositivitySignScan scan = scan_copositivity_signs(A);
    CopositivityChecker checker(A.rows());
    EXPECT_TRUE(checker.are_negative_components_strictly_copositive(A, scan.negative_neighbors));
}

TEST(CopositivityTest, DisconnectedNegativeComponentsPassSeparately) {
    matrix_int A(6, 6);
    for (size_t row = 0; row < 6; ++row) {
        for (size_t column = 0; column < 6; ++column) A(row, column) = integer(row == column ? 4 : 7);
    }
    A(0, 2) = A(2, 0) = integer(-1);
    A(2, 4) = A(4, 2) = integer(-1);
    A(1, 3) = A(3, 1) = integer(-1);
    A(3, 5) = A(5, 3) = integer(-1);

    const CopositivitySignScan scan = scan_copositivity_signs(A);
    CopositivityChecker checker(A.rows());
    EXPECT_TRUE(checker.are_negative_components_strictly_copositive(A, scan.negative_neighbors));
}

TEST(CopositivityTest, DisconnectedNegativeComponentRejectsWholeMatrix) {
    matrix_int A(5, 5);
    for (size_t row = 0; row < 5; ++row) {
        for (size_t column = 0; column < 5; ++column) A(row, column) = integer(row == column ? 1 : 9);
    }
    A(0, 3) = A(3, 0) = integer(-2);
    A(1, 4) = A(4, 1) = integer(-1);

    const CopositivitySignScan scan = scan_copositivity_signs(A);
    CopositivityChecker checker(A.rows());
    EXPECT_FALSE(checker.are_negative_components_strictly_copositive(A, scan.negative_neighbors));
}

TEST(CopositivityTest, NegativeComponentCanContainHighestSupportedIndex) {
    matrix_int A(64, 64);
    for (size_t row = 0; row < 64; ++row) {
        for (size_t column = 0; column < 64; ++column) A(row, column) = integer(row == column ? 1 : 2);
    }
    A(0, 63) = A(63, 0) = integer(-2);

    const CopositivitySignScan scan = scan_copositivity_signs(A);
    EXPECT_EQ(scan.negative_neighbors[0], bs64::single_bit_at_pos(63));
    EXPECT_EQ(scan.negative_neighbors[63], bs64::single_bit_at_pos(0));
    CopositivityChecker checker(A.rows());
    EXPECT_FALSE(checker.are_negative_components_strictly_copositive(A, scan.negative_neighbors));
}

TEST(CopositivityTest, HadelerEnumerationIncludesBit63) {
    matrix_int A;
    A.set_identity(64);
    A(63, 63) = integer(0);

    CopositivityChecker checker(A.rows());
    EXPECT_FALSE(checker.is_strictly_copositive(A));
}

TEST(CopositivityTest, MultiwordSignScanCrossesTwoWordBoundaries) {
    constexpr size_t dimension = 129;
    matrix_int A(dimension, dimension);
    for (size_t i = 0; i < dimension; ++i) A(i, i) = integer(3);
    A(63, 64) = A(64, 63) = integer(-1);
    A(64, 128) = A(128, 64) = integer(-2);

    const CopositivitySignScanMultiword scan = scan_copositivity_signs_multiword(A);

    EXPECT_TRUE(scan.all_diagonal_positive);
    EXPECT_TRUE(scan.has_negative_off_diagonal);
    EXPECT_TRUE(scan.negative_neighbors[63].is_set_at_pos(64));
    EXPECT_TRUE(scan.negative_neighbors[64].is_set_at_pos(63));
    EXPECT_TRUE(scan.negative_neighbors[64].is_set_at_pos(128));
    EXPECT_TRUE(scan.negative_neighbors[128].is_set_at_pos(64));
    EXPECT_FALSE(scan.all_negative_part_row_sums_positive);
}

TEST(CopositivityTest, MultiwordHadelerRejectsBadPairAcrossBit63Boundary) {
    matrix_int A;
    A.set_identity(65);
    A(63, 64) = A(64, 63) = integer(-2);

    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, MultiwordHadelerReachesThirdWordAtCardinalityThree) {
    matrix_int A;
    A.set_identity(129);
    for (const size_t position : {0u, 64u, 128u}) A(position, position) = integer(2);
    A(0, 64) = A(64, 0) = integer(-1);
    A(0, 128) = A(128, 0) = integer(-1);
    A(64, 128) = A(128, 64) = integer(-1);

    EXPECT_FALSE(is_strictly_copositive(A));
}

TEST(CopositivityTest, MultiwordNegativeComponentsCrossWordBoundaries) {
    constexpr size_t dimension = 129;
    matrix_int A(dimension, dimension);
    for (size_t row = 0; row < dimension; ++row) {
        for (size_t column = 0; column < dimension; ++column) A(row, column) = integer(row == column ? 2 : 5);
    }
    A(63, 64) = A(64, 63) = integer(-1);
    A(64, 128) = A(128, 64) = integer(-1);
    A(63, 128) = A(128, 63) = integer(0);

    const CopositivitySignScanMultiword scan = scan_copositivity_signs_multiword(A);
    CopositivityChecker checker(dimension);

    EXPECT_TRUE(checker.are_negative_components_strictly_copositive(A, scan.negative_neighbors));
}

TEST(CopositivityTest, MultiwordConnectedNegativeGraphUsesWholeHadelerMatrix) {
    constexpr size_t dimension = 65;
    matrix_int A(dimension, dimension);
    for (size_t i = 0; i < dimension; ++i) A(i, i) = integer(1);
    for (size_t i = 0; i + 1 < dimension; ++i) A(i, i + 1) = A(i + 1, i) = integer(-1);

    const CopositivitySignScanMultiword scan = scan_copositivity_signs_multiword(A);
    CopositivityChecker checker(dimension);

    EXPECT_FALSE(checker.are_negative_components_strictly_copositive(A, scan.negative_neighbors));
}

TEST(CopositivityTest, DirectSmallCriteriaCoverBoundaryCases) {
    EXPECT_TRUE(is_strictly_copositive_1x1(integer(1)));
    EXPECT_FALSE(is_strictly_copositive_1x1(integer(0)));

    EXPECT_TRUE(is_strictly_copositive_2x2(integer(1), integer(2), integer(1)));
    EXPECT_TRUE(is_strictly_copositive_2x2(integer(2), integer(-1), integer(2)));
    EXPECT_FALSE(is_strictly_copositive_2x2(integer(1), integer(-1), integer(1)));
    EXPECT_FALSE(is_strictly_copositive_2x2(integer(1), integer(-2), integer(1)));

    EXPECT_TRUE(is_strictly_copositive_3x3(
        integer(1), integer(1), integer(1), integer(1), integer(1), integer(1)));
    EXPECT_TRUE(is_strictly_copositive_3x3(
        integer(1), integer(2), integer(2), integer(1), integer(2), integer(1)));
    EXPECT_FALSE(is_strictly_copositive_3x3(
        integer(2), integer(-1), integer(-1), integer(2), integer(-1), integer(2)));
    EXPECT_FALSE(is_strictly_copositive_3x3(
        integer(1), integer(-2), integer(0), integer(1), integer(0), integer(1)));
}
