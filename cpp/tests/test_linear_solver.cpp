#include <gtest/gtest.h>
#include <rational_linalg/linear_solver.hpp>
#include <cmath>

using namespace rational_linalg;

// Test BareissGauss solver
TEST(LinearSolverTest, BareissGauss2x2) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 3, 2x - y = 0
    // Solution: x = 1, y = 2
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 3;
    Ab(1, 0) = 2; Ab(1, 1) = -1; Ab(1, 2) = 0;
    
    BareissGauss<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_TRUE(success);
    EXPECT_EQ(x.rows(), 2);
    EXPECT_EQ(x.cols(), 1);
    EXPECT_EQ(x(0, 0), 1);
    EXPECT_EQ(x(1, 0), 2);
}

TEST(LinearSolverTest, BareissGauss3x3) {
    Matrix<fraction> Ab(3, 4);
    // System: x + y + z = 6, 2x + y - z = 1, x - y + z = 2
    // Solution: x = 1, y = 2, z = 3
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 1; Ab(0, 3) = 6;
    Ab(1, 0) = 2; Ab(1, 1) = 1; Ab(1, 2) = -1; Ab(1, 3) = 1;
    Ab(2, 0) = 1; Ab(2, 1) = -1; Ab(2, 2) = 1; Ab(2, 3) = 2;
    
    BareissGauss<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_TRUE(success);
    EXPECT_EQ(x.rows(), 3);
    EXPECT_EQ(x.cols(), 1);
    EXPECT_EQ(x(0, 0), 1);
    EXPECT_EQ(x(1, 0), 2);
    EXPECT_EQ(x(2, 0), 3);
}

TEST(LinearSolverTest, BareissGaussSingular) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 3, 2x + 2y = 6 (linearly dependent)
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 3;
    Ab(1, 0) = 2; Ab(1, 1) = 2; Ab(1, 2) = 6;
    
    BareissGauss<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_FALSE(success); // Singular system
}

TEST(LinearSolverTest, BareissGaussNegativeSolution) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 1, x - y = 3
    // Solution: x = 2, y = -1 (negative solution should fail)
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 1;
    Ab(1, 0) = 1; Ab(1, 1) = -1; Ab(1, 2) = 3;
    
    BareissGauss<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_FALSE(success); // Negative solution not allowed
}

// Test GaussDouble solver
TEST(LinearSolverTest, GaussDouble2x2) {
    Matrix<double> Ab(2, 3);
    // System: x + y = 3, 2x - y = 0
    // Solution: x = 1, y = 2
    Ab(0, 0) = 1.0; Ab(0, 1) = 1.0; Ab(0, 2) = 3.0;
    Ab(1, 0) = 2.0; Ab(1, 1) = -1.0; Ab(1, 2) = 0.0;
    
    GaussDouble solver(Ab);
    Matrix<double> x;
    bool success = solver.solve(x);
    
    EXPECT_TRUE(success);
    EXPECT_EQ(x.rows(), 2);
    EXPECT_EQ(x.cols(), 1);
    EXPECT_NEAR(x(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(x(1, 0), 2.0, 1e-10);
}

TEST(LinearSolverTest, GaussDouble3x3) {
    Matrix<double> Ab(3, 4);
    // System: x + y + z = 6, 2x + y - z = 1, x - y + z = 2
    // Solution: x = 1, y = 2, z = 3
    Ab(0, 0) = 1.0; Ab(0, 1) = 1.0; Ab(0, 2) = 1.0; Ab(0, 3) = 6.0;
    Ab(1, 0) = 2.0; Ab(1, 1) = 1.0; Ab(1, 2) = -1.0; Ab(1, 3) = 1.0;
    Ab(2, 0) = 1.0; Ab(2, 1) = -1.0; Ab(2, 2) = 1.0; Ab(2, 3) = 2.0;
    
    GaussDouble solver(Ab);
    Matrix<double> x;
    bool success = solver.solve(x);
    
    EXPECT_TRUE(success);
    EXPECT_EQ(x.rows(), 3);
    EXPECT_EQ(x.cols(), 1);
    EXPECT_NEAR(x(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(x(1, 0), 2.0, 1e-10);
    EXPECT_NEAR(x(2, 0), 3.0, 1e-10);
}

TEST(LinearSolverTest, GaussDoubleSingular) {
    Matrix<double> Ab(2, 3);
    // System: x + y = 3, 2x + 2y = 6 (linearly dependent)
    Ab(0, 0) = 1.0; Ab(0, 1) = 1.0; Ab(0, 2) = 3.0;
    Ab(1, 0) = 2.0; Ab(1, 1) = 2.0; Ab(1, 2) = 6.0;
    
    GaussDouble solver(Ab);
    Matrix<double> x;
    bool success = solver.solve(x);
    
    EXPECT_FALSE(success); // Singular system
}

// Test GaussRational solver
TEST(LinearSolverTest, GaussRational2x2) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 3, 2x - y = 0
    // Solution: x = 1, y = 2
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 3;
    Ab(1, 0) = 2; Ab(1, 1) = -1; Ab(1, 2) = 0;
    
    GaussRational<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_TRUE(success);
    EXPECT_EQ(x.rows(), 2);
    EXPECT_EQ(x.cols(), 1);
    EXPECT_EQ(x(0, 0), 1);
    EXPECT_EQ(x(1, 0), 2);
}

TEST(LinearSolverTest, GaussRational3x3) {
    Matrix<fraction> Ab(3, 4);
    // System: x + y + z = 6, 2x + y - z = 1, x - y + z = 2
    // Solution: x = 1, y = 2, z = 3
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 1; Ab(0, 3) = 6;
    Ab(1, 0) = 2; Ab(1, 1) = 1; Ab(1, 2) = -1; Ab(1, 3) = 1;
    Ab(2, 0) = 1; Ab(2, 1) = -1; Ab(2, 2) = 1; Ab(2, 3) = 2;
    
    GaussRational<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_TRUE(success);
    EXPECT_EQ(x.rows(), 3);
    EXPECT_EQ(x.cols(), 1);
    EXPECT_EQ(x(0, 0), 1);
    EXPECT_EQ(x(1, 0), 2);
    EXPECT_EQ(x(2, 0), 3);
}

// Test solver consistency
TEST(LinearSolverTest, SolverConsistency) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 3, 2x - y = 0
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 3;
    Ab(1, 0) = 2; Ab(1, 1) = -1; Ab(1, 2) = 0;
    
    BareissGauss<fraction> bareiss_solver(Ab);
    Matrix<fraction> x_bareiss;
    bool bareiss_success = bareiss_solver.solve(x_bareiss);
    
    GaussRational<fraction> gauss_solver(Ab);
    Matrix<fraction> x_gauss;
    bool gauss_success = gauss_solver.solve(x_gauss);
    
    EXPECT_EQ(bareiss_success, gauss_success);
    if (bareiss_success && gauss_success) {
        EXPECT_EQ(x_bareiss(0, 0), x_gauss(0, 0));
        EXPECT_EQ(x_bareiss(1, 0), x_gauss(1, 0));
    }
}

// Test with fraction type (GMP)
TEST(LinearSolverTest, BareissGaussRationalGMP) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 3, 2x - y = 0
    Ab(0, 0) = fraction(1); Ab(0, 1) = fraction(1); Ab(0, 2) = fraction(3);
    Ab(1, 0) = fraction(2); Ab(1, 1) = fraction(-1); Ab(1, 2) = fraction(0);
    
    BareissGauss<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_TRUE(success);
    EXPECT_EQ(x.rows(), 2);
    EXPECT_EQ(x.cols(), 1);
    EXPECT_EQ(x(0, 0), fraction(1));
    EXPECT_EQ(x(1, 0), fraction(2));
}

// Test zero solution (should fail)
TEST(LinearSolverTest, ZeroSolutionFails) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 0, 2x - y = 0
    // Solution: x = 0, y = 0 (should fail)
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 0;
    Ab(1, 0) = 2; Ab(1, 1) = -1; Ab(1, 2) = 0;
    
    BareissGauss<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_FALSE(success); // Zero solution not allowed
}

// Test positive solution requirement
TEST(LinearSolverTest, PositiveSolutionRequired) {
    Matrix<fraction> Ab(2, 3);
    // System: x + y = 1, x - y = 3
    // Solution: x = 2, y = -1 (negative y should fail)
    Ab(0, 0) = 1; Ab(0, 1) = 1; Ab(0, 2) = 1;
    Ab(1, 0) = 1; Ab(1, 1) = -1; Ab(1, 2) = 3;
    
    BareissGauss<fraction> solver(Ab);
    Matrix<fraction> x;
    bool success = solver.solve(x);
    
    EXPECT_FALSE(success); // Negative solution not allowed
}

