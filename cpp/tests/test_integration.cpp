#include <gtest/gtest.h>
#include <fracessa/fracessa.hpp>
#include <linalg/matrix_fraction.hpp>

using namespace linalg;

TEST(IntegrationTest, SimpleGame) {
    // 3x3 Rock-Paper-Scissors matrix
    std::vector<fraction> vals = {fraction::zero(), fraction::neg_one(), fraction::one(), fraction::zero(), fraction::neg_one(), fraction::zero()};
    matrix_frc A = create_symmetric(3, vals);
    // Actually symmetric version of RPS is not standard, let's just test a simple 2x2
    
    matrix_frc B = matrix_frc(2, 2);
    B(0, 0) = fraction::zero(); B(0, 1) = fraction::one();
    B(1, 0) = fraction::one(); B(1, 1) = fraction::zero();
    
    fracessa analyzer(B, false, true, true, false, false);
    EXPECT_EQ(analyzer.ess_count_, 1);
}
