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
    
    fracessa analyzer(search_method::safe, B, false, true, false, false);
    EXPECT_EQ(analyzer.ess_count_, 1);
}

TEST(IntegrationTest, FullSupportModeReturnsStableFullSupport) {
    matrix_frc B(2, 2);
    B(0, 0) = fraction::zero(); B(0, 1) = fraction::one();
    B(1, 0) = fraction::one();  B(1, 1) = fraction::zero();

    fracessa analyzer(search_method::safe, B, false, true, true, false);

    ASSERT_EQ(analyzer.ess_count_, 1u);
    ASSERT_EQ(analyzer.candidates_.size(), 1u);
    EXPECT_EQ(analyzer.candidates_[0].support_size, 2u);
    EXPECT_FALSE(analyzer.candidates_[0].multiplier.has_value());
}

TEST(IntegrationTest, CircularFullSupportHasMultiplierOne) {
    matrix_frc B(2, 2);
    B(0, 0) = fraction::zero(); B(0, 1) = fraction::one();
    B(1, 0) = fraction::one();  B(1, 1) = fraction::zero();

    fracessa analyzer(search_method::safe, B, true, true, true, false);

    ASSERT_EQ(analyzer.ess_count_, 1u);
    ASSERT_EQ(analyzer.candidates_.size(), 1u);
    ASSERT_TRUE(analyzer.candidates_[0].multiplier.has_value());
    EXPECT_EQ(*analyzer.candidates_[0].multiplier, 1u);
}

TEST(IntegrationTest, FullSupportModeFallsBackWhenFullSupportIsNotStable) {
    matrix_frc B(2, 2);
    B(0, 0) = fraction::one();  B(0, 1) = fraction::zero();
    B(1, 0) = fraction::zero(); B(1, 1) = fraction::one();

    fracessa analyzer(search_method::safe, B, false, true, true, false);

    ASSERT_EQ(analyzer.ess_count_, 2u);
    ASSERT_EQ(analyzer.candidates_.size(), 3u);
    EXPECT_FALSE(analyzer.candidates_[0].is_ess);
    EXPECT_EQ(analyzer.candidates_[0].support_size, 2u);
}

TEST(IntegrationTest, CircularSymmetricStoresOneRepresentative) {
    std::vector<fraction> half_row = {fraction::one(), fraction(3)};
    matrix_frc B = create_circular_symmetric(5, half_row);

    fracessa analyzer(search_method::safe, B, true, true, false, false);

    ASSERT_EQ(analyzer.ess_count_, 5u);
    ASSERT_EQ(analyzer.candidates_.size(), 1u);

    const auto& candidate = analyzer.candidates_[0];
    EXPECT_EQ(candidate.candidate_id, 1u);
    ASSERT_TRUE(candidate.multiplier.has_value());
    EXPECT_EQ(*candidate.multiplier, 5u);
    EXPECT_TRUE(candidate.is_ess);
    EXPECT_EQ(candidate.support_size, 3u);
}

TEST(IntegrationTest, AutomaticCyclicSymmetryPreservesRepresentedResults) {
    const matrix_frc matrix = create_circular_symmetric(8, {fraction(1), fraction(-3), fraction(1), fraction(-3)});

    fracessa exhaustive(search_method::safe, matrix, false, true, false, false);
    fracessa circular(search_method::safe, matrix, true, true, false, false);

    auto represented_candidates = [](const fracessa& analyzer) {
        size_t count = 0;
        for (const candidate& row : analyzer.candidates_)
            count += row.multiplier.value_or(1);
        return count;
    };

    EXPECT_EQ(circular.ess_count_, exhaustive.ess_count_);
    EXPECT_EQ(represented_candidates(circular), represented_candidates(exhaustive));
    ASSERT_EQ(circular.candidates_.size(), 2u);
    EXPECT_LT(circular.candidates_.size(), exhaustive.candidates_.size());

    EXPECT_EQ(circular.candidates_[0].support, 0b00000011u);
    EXPECT_EQ(circular.candidates_[1].support, 0b00001001u);
    for (const candidate& row : circular.candidates_) {
        ASSERT_TRUE(row.multiplier.has_value());
        EXPECT_EQ(*row.multiplier, 8u);
    }
    EXPECT_EQ(circular.candidates_[1].vector(0, 0), fraction(1, 2));
    EXPECT_EQ(circular.candidates_[1].vector(3, 0), fraction(1, 2));
}
