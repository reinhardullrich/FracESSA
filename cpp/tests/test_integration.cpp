#include <gtest/gtest.h>
#include <fracessa/fracessa.hpp>
#include <linalg/matrix_fraction.hpp>

using namespace linalg;

/*
 * End-to-end smoke test:
 * build analyzer from a tiny game matrix and verify ESS count output.
 */

TEST(IntegrationTest, SimpleGame) {
    // Test a simple 2x2 matrix.
    matrix_frc B = matrix_frc(2, 2);
    B(0, 0) = fraction::zero(); B(0, 1) = fraction::one();
    B(1, 0) = fraction::one(); B(1, 1) = fraction::zero();
    
    fracessa analyzer(B, false, true, true, false, false);
    EXPECT_EQ(analyzer.ess_count_, 1);
}

TEST(IntegrationTest, FullSupportModeReturnsStableFullSupport) {
    matrix_frc B(2, 2);
    B(0, 0) = fraction::zero(); B(0, 1) = fraction::one();
    B(1, 0) = fraction::one();  B(1, 1) = fraction::zero();

    fracessa analyzer(B, false, true, true, true, false);

    ASSERT_EQ(analyzer.ess_count_, 1u);
    ASSERT_EQ(analyzer.candidates_.size(), 1u);
    EXPECT_EQ(analyzer.candidates_[0].support_size, 2u);
}

TEST(IntegrationTest, FullSupportModeFallsBackWhenFullSupportIsNotStable) {
    matrix_frc B(2, 2);
    B(0, 0) = fraction::one();  B(0, 1) = fraction::zero();
    B(1, 0) = fraction::zero(); B(1, 1) = fraction::one();

    fracessa analyzer(B, false, true, true, true, false);

    ASSERT_EQ(analyzer.ess_count_, 2u);
    ASSERT_EQ(analyzer.candidates_.size(), 3u);
    EXPECT_FALSE(analyzer.candidates_[0].is_ess);
    EXPECT_EQ(analyzer.candidates_[0].support_size, 2u);
}

TEST(IntegrationTest, CircularSymmetricOrbitShiftReference) {
    std::vector<fraction> half_row = {fraction::one(), fraction(3)};
    matrix_frc B = create_circular_symmetric(5, half_row);

    fracessa analyzer(B, true, true, true, false, false);

    ASSERT_EQ(analyzer.ess_count_, 5u);
    ASSERT_EQ(analyzer.candidates_.size(), 5u);

    for (size_t i = 0; i < analyzer.candidates_.size(); ++i) {
        const auto& c = analyzer.candidates_[i];
        EXPECT_EQ(c.candidate_id, i + 1);
        EXPECT_EQ(c.shift_reference, 1u);
        EXPECT_TRUE(c.is_ess);
        EXPECT_EQ(c.support_size, 3u);
    }
}
