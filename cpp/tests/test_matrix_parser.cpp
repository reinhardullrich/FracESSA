#include <gtest/gtest.h>

#include <fracessa/matrix_parser.hpp>

#include <sstream>
#include <string>

using namespace linalg;

namespace {

std::string join_ones(size_t count) {
    std::ostringstream ss;
    for (size_t i = 0; i < count; ++i) {
        if (i > 0) {
            ss << ",";
        }
        ss << "1";
    }
    return ss.str();
}

} // namespace

TEST(MatrixParserSafeTest, ParsesUpperTriangularSymmetricInput) {
    matrix_frc A;
    bool is_cs = true;

    ASSERT_TRUE(matrix_parser::parse_matrix_string("3#1,2,3,4,5,6", A, is_cs));
    EXPECT_FALSE(is_cs);
    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_EQ(A(0, 0), fraction(1));
    EXPECT_EQ(A(0, 1), fraction(2));
    EXPECT_EQ(A(1, 0), fraction(2));
    EXPECT_EQ(A(2, 2), fraction(6));
}

TEST(MatrixParserSafeTest, ParsesCircularSymmetricInput) {
    matrix_frc A;
    bool is_cs = false;

    ASSERT_TRUE(matrix_parser::parse_matrix_string("5#1,3", A, is_cs));
    EXPECT_TRUE(is_cs);
    EXPECT_EQ(A.rows(), 5);
    EXPECT_EQ(A.cols(), 5);
    EXPECT_EQ(A(0, 0), fraction::zero());
    EXPECT_EQ(A(0, 1), fraction::one());
    EXPECT_EQ(A(0, 2), fraction(3));
    EXPECT_EQ(A(0, 3), fraction(3));
    EXPECT_EQ(A(0, 4), fraction::one());
}

TEST(MatrixParserSafeTest, RejectsMissingHash) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_FALSE(matrix_parser::parse_matrix_string("3,1,2,3,4,5,6", A, is_cs));
}

TEST(MatrixParserSafeTest, RejectsMultipleHashes) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_FALSE(matrix_parser::parse_matrix_string("3#1,2#3", A, is_cs));
}

TEST(MatrixParserSafeTest, RejectsInvalidDimensionRange) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_FALSE(matrix_parser::parse_matrix_string("0#1", A, is_cs));
    EXPECT_FALSE(matrix_parser::parse_matrix_string("64#1", A, is_cs));
}

TEST(MatrixParserSafeTest, RejectsUnexpectedValueCount) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_FALSE(matrix_parser::parse_matrix_string("4#1,2,3", A, is_cs));
}

TEST(MatrixParserUnsafeTest, BypassesSafeDimensionGuard) {
    matrix_frc A;
    bool is_cs = false;

    // n=64 would be rejected by safe parser; unsafe parser accepts and builds CS matrix.
    const std::string payload = "64#" + join_ones(32);
    matrix_parser::parse_matrix_string_unsafe(payload, A, is_cs);

    EXPECT_TRUE(is_cs);
    EXPECT_EQ(A.rows(), 64);
    EXPECT_EQ(A.cols(), 64);
    EXPECT_EQ(A(0, 0), fraction::zero());
    EXPECT_EQ(A(0, 1), fraction::one());
}

TEST(MatrixParserUnsafeTest, FallsBackToSymmetricWhenCountNotCs) {
    matrix_frc A;
    bool is_cs = true;

    matrix_parser::parse_matrix_string_unsafe("3#1,2,3,4,5,6", A, is_cs);

    EXPECT_FALSE(is_cs);
    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A(0, 0), fraction::one());
    EXPECT_EQ(A(1, 0), fraction::two());
}

