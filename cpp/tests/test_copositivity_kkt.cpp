#include <gtest/gtest.h>
#include <linalg/copositive_integer_kkt.hpp>

#include <initializer_list>
#include <random>
#include <stdexcept>

using namespace linalg;

namespace {

matrix_int symmetric_matrix(size_t dimension, std::initializer_list<slong> upper_triangle)
{
    matrix_int matrix(dimension, dimension);
    auto value = upper_triangle.begin();
    for (size_t row = 0; row < dimension; ++row) {
        for (size_t column = row; column < dimension; ++column) {
            matrix(row, column) = integer(*value++);
            matrix(column, row) = matrix(row, column);
        }
    }
    return matrix;
}

TEST(CopositivityKktTest, ClassifiesInteriorAndBoundaryMinimaExactly)
{
    EXPECT_TRUE(is_strictly_copositive_kkt(symmetric_matrix(1, {1})));
    EXPECT_FALSE(is_strictly_copositive_kkt(symmetric_matrix(1, {0})));
    EXPECT_TRUE(is_strictly_copositive_kkt(symmetric_matrix(2, {1, 0, 1})));
    EXPECT_FALSE(is_strictly_copositive_kkt(symmetric_matrix(2, {1, -1, 1})));
}

TEST(CopositivityKktTest, UsesFirstOrderCandidatesWithoutSecondOrderFiltering)
{
    matrix_int identity;
    identity.set_identity(4);
    EXPECT_TRUE(is_strictly_copositive_kkt(identity));

    matrix_int negative_value(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column)
            negative_value(row, column) = integer(row == column ? 5 : -2);
    }
    EXPECT_FALSE(is_strictly_copositive_kkt(negative_value));

    matrix_int zero_value(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column)
            zero_value(row, column) = integer(row == column ? 3 : -1);
    }
    EXPECT_FALSE(is_strictly_copositive_kkt(zero_value));
}

TEST(CopositivityKktTest, PreservesArbitraryPrecisionScaling)
{
    integer scale;
    ASSERT_EQ(fmpz_set_str(scale.native_handle(), "123456789012345678901234567890", 10), 0);

    matrix_int matrix(4, 4);
    for (size_t row = 0; row < 4; ++row) {
        for (size_t column = 0; column < 4; ++column) {
            matrix(row, column) = scale;
            matrix(row, column).multiply(row == column ? 3 : 1);
            if (row != column) matrix(row, column).negate();
        }
    }
    EXPECT_FALSE(is_strictly_copositive_kkt(matrix));
}

TEST(CopositivityKktTest, UsesTheExistingMultiwordSupportHarness)
{
    matrix_int strictly_copositive(65, 65);
    for (size_t row = 0; row < 65; ++row) {
        for (size_t column = 0; column < 65; ++column) strictly_copositive(row, column) = integer(1);
    }
    EXPECT_TRUE(is_strictly_copositive_kkt(strictly_copositive));

    // Without immediate rejection at the final singleton, the remaining identity supports would be infeasible to enumerate.
    matrix_int not_strictly_copositive;
    not_strictly_copositive.set_identity(65);
    not_strictly_copositive(64, 64) = integer(0);
    EXPECT_FALSE(is_strictly_copositive_kkt(not_strictly_copositive));
}

TEST(CopositivityKktTest, RejectsInvalidMatrixShapes)
{
    matrix_int empty;
    EXPECT_THROW(is_strictly_copositive_kkt(empty), std::invalid_argument);

    matrix_int non_square(2, 3);
    EXPECT_THROW(is_strictly_copositive_kkt(non_square), std::invalid_argument);

    matrix_int asymmetric;
    asymmetric.set_identity(2);
    asymmetric(0, 1) = integer(1);
    EXPECT_THROW(is_strictly_copositive_kkt(asymmetric), std::invalid_argument);
}

TEST(CopositivityKktTest, MatchesHadelerOnDeterministicSmallCorpus)
{
    std::mt19937 random(0xF12A55A);
    std::uniform_int_distribution<int> diagonal(1, 6);
    std::uniform_int_distribution<int> off_diagonal(-5, 5);

    for (size_t dimension = 4; dimension <= 7; ++dimension) {
        for (size_t sample = 0; sample < 100; ++sample) {
            matrix_int matrix(dimension, dimension);
            for (size_t row = 0; row < dimension; ++row) {
                for (size_t column = row; column < dimension; ++column) {
                    matrix(row, column) = integer(row == column ? diagonal(random) : off_diagonal(random));
                    matrix(column, row) = matrix(row, column);
                }
            }

            EXPECT_EQ(is_strictly_copositive_kkt(matrix), is_strictly_copositive(matrix))
                << "dimension=" << dimension << ", sample=" << sample << '\n' << matrix.to_pretty_string();
        }
    }
}

} // namespace
