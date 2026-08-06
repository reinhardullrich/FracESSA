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
 * These keep basic storage, indexing, and factory invariants stable.
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

TEST(FindCandidateFastTest, UsesExactPrecisionSpanCutoff) {
    const auto fallback = [](const fraction& diagonal_0, const fraction& off_diagonal, const fraction& diagonal_1) {
        matrix_frc A(2, 2);
        A(0, 0) = diagonal_0;   A(0, 1) = off_diagonal;
        A(1, 0) = off_diagonal; A(1, 1) = diagonal_1;
        find_candidate_safe safe(A);
        find_candidate_fast fast(A);
        fast.convert_game_matrix(safe);
        return fast.safe_fallback_reason();
    };

    EXPECT_EQ(fallback(fraction::zero(), fraction(999'999'999), fraction::zero()), safe_fallback::none);
    EXPECT_EQ(fallback(fraction::zero(), fraction(1'000'000'000), fraction::zero()), safe_fallback::none);
    EXPECT_EQ(fallback(fraction::zero(), fraction("1/1000000000"), fraction::zero()), safe_fallback::none);
    EXPECT_EQ(fallback(fraction::one(), fraction(1'000'000'000), fraction::one()), safe_fallback::precision_span);
    EXPECT_EQ(fallback(fraction(1'000'000'000), fraction(1'000'000'001), fraction(1'000'000'000)),
              safe_fallback::precision_span);
}

TEST(FindCandidateTest, UsesIntegerPrecisionSpanAfterRemovingCommonScale) {
    const auto fallback = [](const fraction& diagonal_0, const fraction& off_diagonal, const fraction& diagonal_1) {
        matrix_frc A(2, 2);
        A(0, 0) = diagonal_0;   A(0, 1) = off_diagonal;
        A(1, 0) = off_diagonal; A(1, 1) = diagonal_1;
        find_candidate_safe safe(A);
        find_candidate_test test(A);
        test.convert_game_matrix(safe);
        return test.safe_fallback_reason();
    };

    EXPECT_EQ(fallback(fraction::zero(), fraction(999'999'999), fraction::zero()), safe_fallback::none);
    EXPECT_EQ(fallback(fraction::zero(), fraction(1'000'000'000), fraction::zero()), safe_fallback::none);
    EXPECT_EQ(fallback(fraction::zero(), fraction("1/1000000000"), fraction::zero()), safe_fallback::none);
    EXPECT_EQ(fallback(fraction::one(), fraction(1'000'000'000), fraction::one()), safe_fallback::precision_span);
    EXPECT_EQ(fallback(fraction(1'000'000'000), fraction(1'000'000'001), fraction(1'000'000'000)),
              safe_fallback::precision_span);
}

TEST(FindCandidateFastAndTest, IgnoreOnlyExactZeroRowsDuringEquilibration) {
    matrix_frc harmless_zero_row(3, 3);
    harmless_zero_row(1, 2) = harmless_zero_row(2, 1) = fraction::one();
    find_candidate_safe harmless_safe(harmless_zero_row);
    find_candidate_fast harmless_fast(harmless_zero_row);
    find_candidate_test harmless_test(harmless_zero_row);
    harmless_fast.convert_game_matrix(harmless_safe);
    harmless_test.convert_game_matrix(harmless_safe);

    EXPECT_EQ(harmless_fast.safe_fallback_reason(), safe_fallback::none);
    EXPECT_EQ(harmless_test.safe_fallback_reason(), safe_fallback::none);

    // Coordinate zero is harmlessly empty, while the remaining K_1,15 adjacency block makes BIN equilibration produce an invalid
    // scale. The zero row must not hide that independent failure.
    matrix_frc zero_row_with_bad_active_block(17, 17);
    for (size_t leaf = 2; leaf < 17; ++leaf) {
        zero_row_with_bad_active_block(1, leaf) = zero_row_with_bad_active_block(leaf, 1) = fraction::one();
    }
    find_candidate_safe bad_safe(zero_row_with_bad_active_block);
    find_candidate_fast bad_fast(zero_row_with_bad_active_block);
    find_candidate_test bad_test(zero_row_with_bad_active_block);
    bad_fast.convert_game_matrix(bad_safe);
    bad_test.convert_game_matrix(bad_safe);

    EXPECT_EQ(bad_fast.safe_fallback_reason(), safe_fallback::equilibration_invalid);
    EXPECT_EQ(bad_test.safe_fallback_reason(), safe_fallback::equilibration_invalid);

    // The Fork Graph remains finite but does not satisfy the BIN convergence test within 100 iterations.
    matrix_frc nonconvergent_game(5, 5);
    nonconvergent_game(0, 4) = nonconvergent_game(4, 0) = fraction::one();
    nonconvergent_game(1, 4) = nonconvergent_game(4, 1) = fraction::one();
    nonconvergent_game(2, 3) = nonconvergent_game(3, 2) = fraction::one();
    nonconvergent_game(3, 4) = nonconvergent_game(4, 3) = fraction::one();
    find_candidate_safe nonconvergent_safe(nonconvergent_game);
    find_candidate_fast nonconvergent_fast(nonconvergent_game);
    find_candidate_test nonconvergent_test(nonconvergent_game);
    nonconvergent_fast.convert_game_matrix(nonconvergent_safe);
    nonconvergent_test.convert_game_matrix(nonconvergent_safe);

    EXPECT_EQ(nonconvergent_fast.safe_fallback_reason(), safe_fallback::equilibration_non_convergence);
    EXPECT_EQ(nonconvergent_test.safe_fallback_reason(), safe_fallback::equilibration_non_convergence);
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

    EXPECT_EQ(fast.safe_fallback_reason(), safe_fallback::none);
    EXPECT_EQ(test.safe_fallback_reason(), safe_fallback::none);
    EXPECT_TRUE(fast.find(bitset64{7}, 3));
    EXPECT_TRUE(test.find(bitset64{7}, 3));
}

TEST(FindCandidateFastAndTest, SolveNonsingularZeroDiagonalTwoByTwoPivot) {
    /*
     * With strategy 0 as reference, the full-support reduced system is
     *
     *     H = [0 1],      r = [-1],
     *         [1 0]           [ 1]
     *
     * so y=(1,-1). A scalar-only LDL^T would stop at the zero diagonal and return true as inconclusive. Bunch-Kaufman must use a
     * 2x2 pivot, solve the system, see the negative probability, and return false.
     */
    matrix_frc game(3, 3);
    game(0, 0) = fraction::zero(); game(0, 1) = game(1, 0) = fraction::one(); game(0, 2) = game(2, 0) = fraction::neg_one();
    game(1, 1) = fraction::two();  game(1, 2) = game(2, 1) = fraction::one(); game(2, 2) = fraction(-2);

    find_candidate_safe safe(game);
    find_candidate_fast fast(game);
    find_candidate_test test(game);
    fast.convert_game_matrix(safe);
    test.convert_game_matrix(safe);

    ASSERT_EQ(fast.safe_fallback_reason(), safe_fallback::none);
    ASSERT_EQ(test.safe_fallback_reason(), safe_fallback::none);
    EXPECT_FALSE(fast.find(bitset64{7}, 3));
    EXPECT_FALSE(test.find(bitset64{7}, 3));
}

TEST(FindCandidateFastAndTest, RemoveCommonDenominatorAndNormalizeGameOnce) {
    /*
     * Before the final common scale, support {0,1,2}, with strategy 0 as reference, has
     *
     *     H = diag(10^-8, 1),     y = (1/4, 1/4).
     *
     * Multiplying the complete game by 10^-50 must not make every reduced pivot fall below the absolute cutoff. Fast and test
     * search remove that common denominator and normalize the integer game once, then the outside payoff still rejects the support.
     */
    const fraction common_scale("1/100000000000000000000000000000000000000000000000000");
    const auto scaled = [&](const char* value) { return fraction(value) * common_scale; };
    matrix_frc game(4, 4);
    game(0, 0) = fraction::zero();
    game(0, 1) = game(1, 0) = scaled("-1/400000000");
    game(0, 2) = game(2, 0) = scaled("-1/4");
    game(0, 3) = game(3, 0) = scaled("-1/10");
    game(1, 1) = scaled("1/200000000");
    game(1, 2) = game(2, 1) = scaled("-100000001/400000000");
    game(1, 3) = game(3, 1) = fraction::zero();
    game(2, 2) = scaled("1/2");
    game(2, 3) = game(3, 2) = fraction::zero();
    game(3, 3) = fraction::zero();

    find_candidate_safe safe(game);
    find_candidate_fast fast(game);
    find_candidate_test test(game);
    fast.convert_game_matrix(safe);
    test.convert_game_matrix(safe);

    ASSERT_EQ(fast.safe_fallback_reason(), safe_fallback::none);
    ASSERT_EQ(test.safe_fallback_reason(), safe_fallback::none);
    EXPECT_FALSE(fast.find(bitset64{7}, 3));
    EXPECT_FALSE(test.find(bitset64{7}, 3));
}
