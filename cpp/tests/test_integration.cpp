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

TEST(IntegrationTest, FullSupportModeHandlesDimension64) {
    constexpr size_t dimension = bs64::kMaxBitsetDimension;
    matrix_frc B(dimension, dimension);
    for (size_t strategy = 0; strategy < dimension; ++strategy) B(strategy, strategy) = fraction(-1);

    fracessa analyzer(search_method::safe, B, false, true, true, false);

    ASSERT_EQ(analyzer.ess_count_, 1u);
    ASSERT_EQ(analyzer.candidates_.size(), 1u);
    EXPECT_EQ(analyzer.candidate_structure_[dimension], 1u);
    EXPECT_EQ(analyzer.ess_structure_[dimension], 1u);
    EXPECT_EQ(analyzer.candidates_[0].support, ~bitset64{0});
    EXPECT_EQ(analyzer.candidates_[0].support_size, dimension);
}

TEST(IntegrationTest, CircularSearchHandlesDimension64) {
    constexpr size_t dimension = bs64::kMaxBitsetDimension;
    const matrix_frc B = create_circular_symmetric(dimension, std::vector<fraction>(dimension / 2, fraction(-1)));

    fracessa analyzer(search_method::safe, B, true, true, false, false);

    EXPECT_EQ(analyzer.candidate_count_, dimension);
    EXPECT_EQ(analyzer.ess_count_, dimension);
    EXPECT_EQ(analyzer.candidate_structure_[1], dimension);
    EXPECT_EQ(analyzer.ess_structure_[1], dimension);
    ASSERT_EQ(analyzer.candidates_.size(), 1u);
    ASSERT_TRUE(analyzer.candidates_[0].multiplier.has_value());
    EXPECT_EQ(*analyzer.candidates_[0].multiplier, dimension);
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

TEST(IntegrationTest, EarlyCopositivityDecisionsNameTheirCertificate) {
    matrix_frc diagonal_witness(2, 2);
    diagonal_witness(0, 0) = fraction::one(); diagonal_witness(0, 1) = fraction::one();
    diagonal_witness(1, 0) = fraction::one(); diagonal_witness(1, 1) = fraction::two();

    fracessa diagonal_analyzer(search_method::safe, diagonal_witness, false, true, false, false);
    ASSERT_FALSE(diagonal_analyzer.candidates_.empty());
    EXPECT_EQ(diagonal_analyzer.candidates_[0].stability, "F_not_copos_small");

    matrix_frc all_ones_witness(4, 4);
    all_ones_witness(0, 0) = fraction(0);     all_ones_witness(0, 1) = fraction(0);
    all_ones_witness(0, 2) = fraction(0);     all_ones_witness(0, 3) = fraction(0);
    all_ones_witness(1, 0) = fraction(0);     all_ones_witness(1, 1) = fraction(-53);
    all_ones_witness(1, 2) = fraction(33, 2); all_ones_witness(1, 3) = fraction(73, 2);
    all_ones_witness(2, 0) = fraction(0);     all_ones_witness(2, 1) = fraction(33, 2);
    all_ones_witness(2, 2) = fraction(-85, 2); all_ones_witness(2, 3) = fraction(26);
    all_ones_witness(3, 0) = fraction(0);     all_ones_witness(3, 1) = fraction(73, 2);
    all_ones_witness(3, 2) = fraction(26);    all_ones_witness(3, 3) = fraction(-125, 2);

    fracessa all_ones_analyzer(search_method::safe, all_ones_witness, false, true, false, false);
    ASSERT_FALSE(all_ones_analyzer.candidates_.empty());
    EXPECT_EQ(all_ones_analyzer.candidates_[0].stability, "F_not_copos_small");

    matrix_frc nonnegative_off_diagonal(3, 3);
    nonnegative_off_diagonal(0, 0) = fraction::one(); nonnegative_off_diagonal(0, 1) = fraction::one();
    nonnegative_off_diagonal(0, 2) = fraction::one(); nonnegative_off_diagonal(1, 0) = fraction::one();
    nonnegative_off_diagonal(1, 1) = fraction::zero(); nonnegative_off_diagonal(1, 2) = fraction::zero();
    nonnegative_off_diagonal(2, 0) = fraction::one(); nonnegative_off_diagonal(2, 1) = fraction::zero();
    nonnegative_off_diagonal(2, 2) = fraction::zero();

    fracessa nonnegative_analyzer(search_method::safe, nonnegative_off_diagonal, false, true, false, false);
    ASSERT_FALSE(nonnegative_analyzer.candidates_.empty());
    EXPECT_EQ(nonnegative_analyzer.candidates_[0].stability, "T_copos_small");

    matrix_frc cone_witness(5, 5);
    cone_witness(1, 1) = fraction(-1);
    cone_witness(1, 2) = fraction(2); cone_witness(2, 1) = fraction(2);
    cone_witness(2, 2) = fraction(-1);
    cone_witness(3, 3) = fraction(-10);
    cone_witness(4, 4) = fraction(-10);

    fracessa cone_analyzer(search_method::safe, cone_witness, false, true, false, false);
    ASSERT_FALSE(cone_analyzer.candidates_.empty());
    EXPECT_EQ(cone_analyzer.candidates_[0].support, 1u);
    EXPECT_EQ(cone_analyzer.candidates_[0].extended_support_size, 5u);
    EXPECT_EQ(cone_analyzer.candidates_[0].stability, "F_not_copos");
}

TEST(IntegrationTest, CircularSymmetricStoresOneRepresentative) {
    std::vector<fraction> half_row = {fraction::one(), fraction(3)};
    matrix_frc B = create_circular_symmetric(5, half_row);

    fracessa analyzer(search_method::safe, B, true, true, false, false);

    ASSERT_EQ(analyzer.ess_count_, 5u);
    EXPECT_EQ(analyzer.candidate_count_, 5u);
    EXPECT_EQ(analyzer.candidate_structure_[3], 5u);
    EXPECT_EQ(analyzer.ess_structure_[3], 5u);
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
