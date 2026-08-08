#include <gtest/gtest.h>

#include <fracessa/matrix_parser.hpp>

#include <sstream>
#include <stdexcept>
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

TEST(MatrixParserTest, ParsesUpperTriangularSymmetricInput) {
    matrix_frc A;
    bool is_cs = true;

    ASSERT_NO_THROW(matrix_parser::parse_matrix_string("3#1,2,3,4,5,6", A, is_cs));
    EXPECT_FALSE(is_cs);
    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_EQ(A(0, 0), fraction(1));
    EXPECT_EQ(A(0, 1), fraction(2));
    EXPECT_EQ(A(1, 0), fraction(2));
    EXPECT_EQ(A(2, 2), fraction(6));
}

TEST(MatrixParserTest, ParsesCircularSymmetricInput) {
    matrix_frc A;
    bool is_cs = false;

    ASSERT_NO_THROW(matrix_parser::parse_matrix_string("5#1,3", A, is_cs));
    EXPECT_TRUE(is_cs);
    EXPECT_EQ(A.rows(), 5);
    EXPECT_EQ(A.cols(), 5);
    EXPECT_EQ(A(0, 0), fraction::zero());
    EXPECT_EQ(A(0, 1), fraction::one());
    EXPECT_EQ(A(0, 2), fraction(3));
    EXPECT_EQ(A(0, 3), fraction(3));
    EXPECT_EQ(A(0, 4), fraction::one());
}

TEST(MatrixParserTest, ParsesValuesBeyondWindowsLongRange) {
    matrix_frc A;
    bool is_cs = true;

    ASSERT_NO_THROW(matrix_parser::parse_matrix_string("2#2147483648/2147483649,0,1", A, is_cs));
    EXPECT_FALSE(is_cs);
    EXPECT_EQ(A(0, 0).to_string(), "2147483648/2147483649");
    EXPECT_EQ(A(0, 1), fraction::zero());
    EXPECT_EQ(A(1, 1), fraction::one());
}

TEST(MatrixParserTest, ParsesFastAndArbitraryPrecisionValues) {
    matrix_frc A;
    bool is_cs = true;

    ASSERT_NO_THROW(matrix_parser::parse_matrix_string(
        "3#999999999999999999,1000000000000000000,1/1000000000000000000,"
        "-100000000000000000000,-1/-2,1",
        A, is_cs));

    EXPECT_FALSE(is_cs);
    EXPECT_EQ(A(0, 0).to_string(), "999999999999999999");
    EXPECT_EQ(A(0, 1).to_string(), "1000000000000000000");
    EXPECT_EQ(A(0, 2).to_string(), "1/1000000000000000000");
    EXPECT_EQ(A(1, 1).to_string(), "-100000000000000000000");
    EXPECT_EQ(A(1, 2).to_string(), "1/2");
}

TEST(MatrixParserTest, RejectsZeroDenominator) {
    matrix_frc A;
    bool is_cs = false;

    EXPECT_THROW(matrix_parser::parse_matrix_string("2#1/0,0,1", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("2#1/-0,0,1", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("2#1/00,0,1", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("2#1/0 ,0,1", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("2#1/ 0,0,1", A, is_cs), std::invalid_argument);
}

TEST(MatrixParserTest, RejectsMissingHash) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_THROW(matrix_parser::parse_matrix_string("3,1,2,3,4,5,6", A, is_cs), std::invalid_argument);
}

TEST(MatrixParserTest, RejectsMultipleHashes) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_THROW(matrix_parser::parse_matrix_string("3#1,2#3", A, is_cs), std::invalid_argument);
}

TEST(MatrixParserTest, RejectsInvalidDimensionRange) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_THROW(matrix_parser::parse_matrix_string("0#1", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("65#1", A, is_cs), std::invalid_argument);
}

TEST(MatrixParserTest, RejectsNonDecimalDimensionText) {
    matrix_frc A;
    bool is_cs = false;

    EXPECT_THROW(matrix_parser::parse_matrix_string("2x#0,1,0", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("x2#0,1,0", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string(" 2#0,1,0", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("+2#0,1,0", A, is_cs), std::invalid_argument);
    EXPECT_THROW(matrix_parser::parse_matrix_string("2 #0,1,0", A, is_cs), std::invalid_argument);
}

TEST(MatrixParserTest, AcceptsMaximumSearchDimension) {
    matrix_frc A;
    bool is_cs = false;

    const std::string payload = "64#" + join_ones(32);
    ASSERT_NO_THROW(matrix_parser::parse_matrix_string(payload, A, is_cs));

    EXPECT_TRUE(is_cs);
    EXPECT_EQ(A.rows(), 64);
    EXPECT_EQ(A.cols(), 64);
}

TEST(MatrixParserTest, RejectsUnexpectedValueCount) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_THROW(matrix_parser::parse_matrix_string("4#1,2,3", A, is_cs), std::invalid_argument);
}

TEST(MatrixParserTest, RejectsTrailingComma) {
    matrix_frc A;
    bool is_cs = false;
    EXPECT_THROW(matrix_parser::parse_matrix_string("2#1,2,3,", A, is_cs), std::invalid_argument);
}
