#include <gtest/gtest.h>
#include <coposit/parsers/matrix_parser.hpp>
#include <fracessa/fast_candidate_filter.hpp>
#include <fracessa/types.hpp>

using namespace fracessa::numeric;
using namespace fracessa::search;
using namespace fracessa::support;

/*
 * Matrix container and factory tests.
 *
 * These keep basic storage, indexing, and factory invariants stable.
 */

TEST(MatrixDoubleTest, BasicOperations) {
    matrix_dbl A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 2.0;
    A(1, 0) = 3.0; A(1, 1) = 4.0;

    EXPECT_EQ(A.rows(), 2);
    EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A(1, 1), 4.0);
}

TEST(MatrixIntegerTest, PreservesExactFractionMatrixValues) {
    const auto parsed = coposit::parsers::matrix_parser::parse("2#1/2,-2/3,7");

    EXPECT_EQ(parsed.matrix.rows(), 2u);
    EXPECT_EQ(parsed.matrix.cols(), 2u);
    EXPECT_EQ(parsed.denominator.compare(integer(6)), 0);
    EXPECT_EQ(parsed.matrix(0, 0).compare(integer(3)), 0);
    EXPECT_EQ(parsed.matrix(0, 1).compare(integer(-4)), 0);
    EXPECT_EQ(parsed.matrix(1, 0).compare(integer(-4)), 0);
    EXPECT_EQ(parsed.matrix(1, 1).compare(integer(42)), 0);
}

TEST(FastCandidateFilterTest, UsesExactPrecisionSpanCutoff) {
    const auto fallback = [](std::string_view values) {
        const auto game = coposit::parsers::matrix_parser::parse(std::string("2#") + std::string(values));
        fast_candidate_filter fast(game.matrix.rows());
        fast.convert_game_matrix(game.matrix);
        return fast.safe_fallback_reason();
    };

    EXPECT_EQ(fallback("0,999999999,0"), safe_fallback::none);
    EXPECT_EQ(fallback("0,1000000000,0"), safe_fallback::none);
    EXPECT_EQ(fallback("0,1/1000000000,0"), safe_fallback::none);
    EXPECT_EQ(fallback("1,1000000000,1"), safe_fallback::precision_span);
    EXPECT_EQ(fallback("1000000000,1000000001,1000000000"), safe_fallback::precision_span);
}

TEST(FastCandidateFilterTest, IgnoreOnlyExactZeroRowsDuringEquilibration) {
    matrix_int harmless_zero_row(3, 3);
    harmless_zero_row(1, 2) = harmless_zero_row(2, 1) = integer(1);
    fast_candidate_filter harmless_fast(harmless_zero_row.rows());
    harmless_fast.convert_game_matrix(harmless_zero_row);

    EXPECT_EQ(harmless_fast.safe_fallback_reason(), safe_fallback::none);

    // Coordinate zero is harmlessly empty, while the remaining K_1,15 adjacency block makes BIN equilibration produce an invalid
    // scale. The zero row must not hide that independent failure.
    matrix_int zero_row_with_bad_active_block(17, 17);
    for (size_t leaf = 2; leaf < 17; ++leaf) {
        zero_row_with_bad_active_block(1, leaf) = zero_row_with_bad_active_block(leaf, 1) = integer(1);
    }
    fast_candidate_filter bad_fast(zero_row_with_bad_active_block.rows());
    bad_fast.convert_game_matrix(zero_row_with_bad_active_block);

    EXPECT_EQ(bad_fast.safe_fallback_reason(), safe_fallback::equilibration_invalid);

    // The Fork Graph remains finite but does not satisfy the BIN convergence test within 100 iterations.
    matrix_int nonconvergent_game(5, 5);
    nonconvergent_game(0, 4) = nonconvergent_game(4, 0) = integer(1);
    nonconvergent_game(1, 4) = nonconvergent_game(4, 1) = integer(1);
    nonconvergent_game(2, 3) = nonconvergent_game(3, 2) = integer(1);
    nonconvergent_game(3, 4) = nonconvergent_game(4, 3) = integer(1);
    fast_candidate_filter nonconvergent_fast(nonconvergent_game.rows());
    nonconvergent_fast.convert_game_matrix(nonconvergent_game);

    EXPECT_EQ(nonconvergent_fast.safe_fallback_reason(), safe_fallback::equilibration_non_convergence);
}

TEST(FastCandidateFilterTest, SendSmallPivotToExactArithmetic) {
    const auto game = coposit::parsers::matrix_parser::parse(
        "3#-3,1,2,-1000000000001/3000000000000,-1999999999999/3000000000000,-4000000000001/3000000000000");
    fast_candidate_filter fast(game.matrix.rows());
    fast.convert_game_matrix(game.matrix);

    EXPECT_EQ(fast.safe_fallback_reason(), safe_fallback::none);
    EXPECT_TRUE(fast.passes(bitset{7}, 3));
}

TEST(FastCandidateFilterTest, SolveNonsingularZeroDiagonalTwoByTwoPivot) {
    /*
     * With strategy 0 as reference, the full-support reduced system is
     *
     *     H = [0 1],      r = [-1],
     *         [1 0]           [ 1]
     *
     * so y=(1,-1). A scalar-only LDL^T would stop at the zero diagonal and return true as inconclusive. Bunch-Kaufman must use a
     * 2x2 pivot, solve the system, see the negative probability, and return false.
     */
    const auto game = coposit::parsers::matrix_parser::parse("3#0,1,-1,2,1,-2");
    fast_candidate_filter fast(game.matrix.rows());
    fast.convert_game_matrix(game.matrix);

    ASSERT_EQ(fast.safe_fallback_reason(), safe_fallback::none);
    EXPECT_FALSE(fast.passes(bitset{7}, 3));
}

TEST(FastCandidateFilterTest, RemoveCommonDenominatorAndNormalizeGameOnce) {
    /*
     * Before the final common scale, support {0,1,2}, with strategy 0 as reference, has
     *
     *     H = diag(10^-8, 1),     y = (1/4, 1/4).
     *
     * Multiplying the complete game by 10^-50 must not make every reduced pivot fall below the absolute cutoff. Fast search removes
     * that common denominator and normalizes the integer game once, then the outside payoff still rejects the support.
     */
    const auto game = coposit::parsers::matrix_parser::parse(
        "4#0,-1/40000000000000000000000000000000000000000000000000000000000,"
        "-1/400000000000000000000000000000000000000000000000000,-1/1000000000000000000000000000000000000000000000000000,"
        "1/20000000000000000000000000000000000000000000000000000000000,"
        "-100000001/40000000000000000000000000000000000000000000000000000000000,0,"
        "1/200000000000000000000000000000000000000000000000000,0,0");
    fast_candidate_filter fast(game.matrix.rows());
    fast.convert_game_matrix(game.matrix);

    ASSERT_EQ(fast.safe_fallback_reason(), safe_fallback::none);
    EXPECT_FALSE(fast.passes(bitset{7}, 3));
}

TEST(FastCandidateFilterTest, MultiwordOutsideTraversalCrossesBit63And64)
{
    constexpr size_t dimension = 65;
    matrix_int game(dimension, dimension);
    for (size_t strategy = 0; strategy < dimension; ++strategy) game(strategy, strategy) = integer(1);
    game(0, 64) = game(64, 0) = integer(2);

    fast_candidate_filter fast(game.rows());
    fast.convert_game_matrix(game);

    ASSERT_EQ(fast.safe_fallback_reason(), safe_fallback::none);

    bitset_multiword support(dimension);
    support.set_bit_at_pos(0);
    EXPECT_FALSE(fast.passes(support, 1));
}
