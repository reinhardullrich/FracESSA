#include <gtest/gtest.h>
#include <fracessa/find_candidate_fast.hpp>
#include <fracessa/find_candidate_safe.hpp>
#include <fracessa/find_candidate_test.hpp>
#include <linalg/integer.hpp>
#include <linalg/matrix_fraction.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_integer.hpp>

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

TEST(MatrixIntegerTest, PreservesExactFractionMatrixValues) {
    matrix_frc rational(2, 2);
    rational(0, 0) = fraction(1, 2);  rational(0, 1) = fraction(-2, 3);
    rational(1, 0) = fraction(5, 6);  rational(1, 1) = fraction(7);

    matrix_int scaled;
    integer denominator;
    scaled.set_from_fraction_matrix(rational, denominator);

    EXPECT_EQ(scaled.rows(), rational.rows());
    EXPECT_EQ(scaled.cols(), rational.cols());
    EXPECT_GT(denominator.sign(), 0);
    for (size_t row = 0; row < rational.rows(); ++row) {
        for (size_t column = 0; column < rational.cols(); ++column) {
            fraction restored;
            restored.set_ratio(scaled(row, column), denominator);
            EXPECT_EQ(restored, rational(row, column));
        }
    }
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

TEST(FindCandidateFastTest, UsesExactPrecisionSpanCutoff) {
    const auto exceeds_cutoff = [](const fraction& diagonal_0, const fraction& off_diagonal, const fraction& diagonal_1) {
        matrix_frc A(2, 2);
        A(0, 0) = diagonal_0;   A(0, 1) = off_diagonal;
        A(1, 0) = off_diagonal; A(1, 1) = diagonal_1;
        find_candidate_safe safe(A);
        find_candidate_fast fast(A);
        fast.convert_game_matrix(safe);
        return fast.requires_safe_fallback();
    };

    EXPECT_FALSE(exceeds_cutoff(fraction::zero(), fraction(999'999'999), fraction::zero()));
    EXPECT_TRUE(exceeds_cutoff(fraction::zero(), fraction(1'000'000'000), fraction::zero()));
    EXPECT_TRUE(exceeds_cutoff(fraction::zero(), fraction("1/1000000000"), fraction::zero()));
    EXPECT_TRUE(exceeds_cutoff(fraction(1'000'000'000), fraction(1'000'000'001), fraction(1'000'000'000)));
}

TEST(FindCandidateTest, UsesExactPrecisionSpanCutoff) {
    const auto exceeds_cutoff = [](const fraction& diagonal_0, const fraction& off_diagonal, const fraction& diagonal_1) {
        matrix_frc A(2, 2);
        A(0, 0) = diagonal_0;   A(0, 1) = off_diagonal;
        A(1, 0) = off_diagonal; A(1, 1) = diagonal_1;
        find_candidate_safe safe(A);
        find_candidate_test test(A);
        test.convert_game_matrix(safe);
        return test.requires_safe_fallback();
    };

    EXPECT_FALSE(exceeds_cutoff(fraction::zero(), fraction(999'999'999), fraction::zero()));
    EXPECT_TRUE(exceeds_cutoff(fraction::zero(), fraction(1'000'000'000), fraction::zero()));
    EXPECT_TRUE(exceeds_cutoff(fraction::zero(), fraction("1/1000000000"), fraction::zero()));
    EXPECT_TRUE(exceeds_cutoff(fraction(1'000'000'000), fraction(1'000'000'001), fraction(1'000'000'000)));
}

TEST(FindCandidateFastAndTest, SendSmallPivotToExactArithmetic) {
    matrix_frc A(3, 3);
    A(0, 0) = fraction(-3);
    A(0, 1) = A(1, 0) = fraction(1);
    A(0, 2) = A(2, 0) = fraction(2);
    A(1, 1) = fraction("-1000000000001/3000000000000");
    A(1, 2) = A(2, 1) = fraction("-1999999999999/3000000000000");
    A(2, 2) = fraction("-4000000000001/3000000000000");

    find_candidate_fast fast(A);
    find_candidate_safe safe(A);
    find_candidate_test test(A);
    fast.convert_game_matrix(safe);
    test.convert_game_matrix(safe);

    EXPECT_FALSE(fast.requires_safe_fallback());
    EXPECT_FALSE(test.requires_safe_fallback());
    EXPECT_TRUE(fast.find(bitset64{7}, 3));
    EXPECT_TRUE(test.find(bitset64{7}, 3));
}
