#include <gtest/gtest.h>
#include <rational_linalg/matrix.hpp>
#include <limits>

using namespace rational_linalg;

// Test Construction
TEST(MatrixTest, DefaultConstruction) {
    Matrix<fraction> m;
    EXPECT_EQ(m.rows(), 0);
    EXPECT_EQ(m.cols(), 0);
}

TEST(MatrixTest, SizedConstruction) {
    Matrix<fraction> m(3, 4);
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 4);
    // All elements should be zero-initialized
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            EXPECT_EQ(m(i, j), fraction(0));
        }
    }
}

TEST(MatrixTest, CopyConstruction) {
    Matrix<fraction> m1(2, 2);
    m1(0, 0) = 1; m1(0, 1) = 2;
    m1(1, 0) = 3; m1(1, 1) = 4;
    
    Matrix<fraction> m2(m1);
    EXPECT_EQ(m2.rows(), 2);
    EXPECT_EQ(m2.cols(), 2);
    EXPECT_EQ(m2(0, 0), 1);
    EXPECT_EQ(m2(0, 1), 2);
    EXPECT_EQ(m2(1, 0), 3);
    EXPECT_EQ(m2(1, 1), 4);
}

TEST(MatrixTest, MoveConstruction) {
    Matrix<fraction> m1(2, 2);
    m1(0, 0) = 1;
    
    Matrix<fraction> m2(std::move(m1));
    EXPECT_EQ(m2.rows(), 2);
    EXPECT_EQ(m2.cols(), 2);
    EXPECT_EQ(m2(0, 0), 1);
    EXPECT_EQ(m1.rows(), 0); // Moved from
}

// Test Accessors
TEST(MatrixTest, ElementAccess) {
    Matrix<fraction> m(2, 2);
    m(0, 0) = 1;
    m(0, 1) = 2;
    m(1, 0) = 3;
    m(1, 1) = 4;
    
    EXPECT_EQ(m(0, 0), 1);
    EXPECT_EQ(m(0, 1), 2);
    EXPECT_EQ(m(1, 0), 3);
    EXPECT_EQ(m(1, 1), 4);
}

TEST(MatrixTest, ConstElementAccess) {
    Matrix<fraction> m(2, 2);
    m(0, 0) = 1;
    const Matrix<fraction>& cm = m;
    EXPECT_EQ(cm(0, 0), 1);
}

TEST(MatrixTest, DataAccess) {
    Matrix<fraction> m(2, 2);
    m(0, 0) = 1;
    EXPECT_EQ(m.data()[0], 1);
}

// Test Arithmetic Operations


TEST(MatrixTest, MatrixMultiplication) {
    Matrix<fraction> a(2, 3);
    a(0, 0) = 1; a(0, 1) = 2; a(0, 2) = 3;
    a(1, 0) = 4; a(1, 1) = 5; a(1, 2) = 6;
    
    Matrix<fraction> b(3, 2);
    b(0, 0) = 7; b(0, 1) = 8;
    b(1, 0) = 9; b(1, 1) = 10;
    b(2, 0) = 11; b(2, 1) = 12;
    
    Matrix<fraction> c = a * b;
    EXPECT_EQ(c.rows(), 2);
    EXPECT_EQ(c.cols(), 2);
    EXPECT_EQ(c(0, 0), 58);  // 1*7 + 2*9 + 3*11
    EXPECT_EQ(c(0, 1), 64);  // 1*8 + 2*10 + 3*12
    EXPECT_EQ(c(1, 0), 139); // 4*7 + 5*9 + 6*11
    EXPECT_EQ(c(1, 1), 154); // 4*8 + 5*10 + 6*12
}

TEST(MatrixTest, ScalarMultiplication) {
    Matrix<fraction> a(2, 2);
    a(0, 0) = 1; a(0, 1) = 2;
    a(1, 0) = 3; a(1, 1) = 4;
    
    Matrix<fraction> b = a * fraction(2);
    EXPECT_EQ(b(0, 0), 2);
    EXPECT_EQ(b(0, 1), 4);
    EXPECT_EQ(b(1, 0), 6);
    EXPECT_EQ(b(1, 1), 8);
}



TEST(MatrixTest, CompoundAssignmentAddition) {
    Matrix<fraction> a(2, 2);
    a(0, 0) = 1; a(0, 1) = 2;
    a(1, 0) = 3; a(1, 1) = 4;
    
    Matrix<fraction> b(2, 2);
    b(0, 0) = 5; b(0, 1) = 6;
    b(1, 0) = 7; b(1, 1) = 8;
    
    a += b;
    EXPECT_EQ(a(0, 0), 6);
    EXPECT_EQ(a(0, 1), 8);
    EXPECT_EQ(a(1, 0), 10);
    EXPECT_EQ(a(1, 1), 12);
}

// Test Matrix Operations
TEST(MatrixTest, Transpose) {
    Matrix<fraction> a(2, 3);
    a(0, 0) = 1; a(0, 1) = 2; a(0, 2) = 3;
    a(1, 0) = 4; a(1, 1) = 5; a(1, 2) = 6;
    
    Matrix<fraction> b = a.transpose();
    EXPECT_EQ(b.rows(), 3);
    EXPECT_EQ(b.cols(), 2);
    EXPECT_EQ(b(0, 0), 1);
    EXPECT_EQ(b(0, 1), 4);
    EXPECT_EQ(b(1, 0), 2);
    EXPECT_EQ(b(1, 1), 5);
    EXPECT_EQ(b(2, 0), 3);
    EXPECT_EQ(b(2, 1), 6);
}






// Test Factory Functions
TEST(MatrixTest, ZeroFactory) {
    Matrix<fraction> z = Matrix<fraction>::zero(3, 4);
    EXPECT_EQ(z.rows(), 3);
    EXPECT_EQ(z.cols(), 4);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            EXPECT_EQ(z(i, j), 0);
        }
    }
}

TEST(MatrixTest, IdentityFactory) {
    Matrix<fraction> I = Matrix<fraction>::identity(3);
    EXPECT_EQ(I.rows(), 3);
    EXPECT_EQ(I.cols(), 3);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            if (i == j) {
                EXPECT_EQ(I(i, j), 1);
            } else {
                EXPECT_EQ(I(i, j), 0);
            }
        }
    }
}


TEST(MatrixTest, ZeroVectorFactory) {
    Matrix<fraction> z = Matrix<fraction>::Zero(4);
    EXPECT_EQ(z.rows(), 4);
    EXPECT_EQ(z.cols(), 1);
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(z(i, 0), 0);
    }
}

// Test Utility Functions
TEST(MatrixTest, SwapRows) {
    Matrix<fraction> a(2, 2);
    a(0, 0) = 1; a(0, 1) = 2;
    a(1, 0) = 3; a(1, 1) = 4;
    
    a.swap_rows(0, 1);
    EXPECT_EQ(a(0, 0), 3);
    EXPECT_EQ(a(0, 1), 4);
    EXPECT_EQ(a(1, 0), 1);
    EXPECT_EQ(a(1, 1), 2);
}





// Test Comparison


// Test Matrix Factory Functions
TEST(MatrixTest, CreateSymmetric) {
    std::vector<fraction> upper = {1, 2, 3, 4, 5, 6};
    Matrix<fraction> m = create_symmetric<fraction>(3, upper);
    
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 3);
    EXPECT_EQ(m(0, 0), 1);
    EXPECT_EQ(m(0, 1), 2);
    EXPECT_EQ(m(0, 2), 3);
    EXPECT_EQ(m(1, 0), 2); // Symmetric
    EXPECT_EQ(m(1, 1), 4);
    EXPECT_EQ(m(1, 2), 5);
    EXPECT_EQ(m(2, 0), 3); // Symmetric
    EXPECT_EQ(m(2, 1), 5); // Symmetric
    EXPECT_EQ(m(2, 2), 6);
}

TEST(MatrixTest, CreateCircularSymmetric) {
    std::vector<fraction> half_row = {1, 2};
    Matrix<fraction> m = create_circular_symmetric<fraction>(4, half_row);
    
    EXPECT_EQ(m.rows(), 4);
    EXPECT_EQ(m.cols(), 4);
    // Check circular property: m(i, j) = m(0, (j-i) mod n)
    EXPECT_EQ(m(0, 0), 0);
    EXPECT_EQ(m(0, 1), 1);
    EXPECT_EQ(m(0, 2), 2);
    EXPECT_EQ(m(0, 3), 1);
}

// Test Conversions
TEST(MatrixTest, ConvertSmallToRational) {
    Matrix<fraction> a(2, 2);
    a(0, 0) = 1; a(0, 1) = 2;
    a(1, 0) = 3; a(1, 1) = 4;
    
    Matrix<fraction> b = a;  // No conversion needed, both use fraction
    EXPECT_EQ(b.rows(), 2);
    EXPECT_EQ(b.cols(), 2);
    EXPECT_EQ(b(0, 0), fraction(1));
    EXPECT_EQ(b(0, 1), fraction(2));
    EXPECT_EQ(b(1, 0), fraction(3));
    EXPECT_EQ(b(1, 1), fraction(4));
}

TEST(MatrixTest, ConvertToDouble) {
    Matrix<fraction> a(2, 2);
    a(0, 0) = 1; a(0, 1) = 2;
    a(1, 0) = 3; a(1, 1) = 4;
    
    Matrix<double> b = a.to_double();
    EXPECT_EQ(b.rows(), 2);
    EXPECT_EQ(b.cols(), 2);
    EXPECT_DOUBLE_EQ(b(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(b(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(b(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(b(1, 1), 4.0);
}

