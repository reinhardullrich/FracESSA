#include <gtest/gtest.h>

#include <fracessa/supports.hpp>

TEST(SupportsTest, InitializeNonCircularCountsByCardinality) {
    Supports supports(4, false);
    supports.initialize();

    EXPECT_EQ(supports.get_supports(1).size(), 4u);
    EXPECT_EQ(supports.get_supports(2).size(), 6u);
    EXPECT_EQ(supports.get_supports(3).size(), 4u);
    EXPECT_EQ(supports.get_supports(4).size(), 1u);
}

TEST(SupportsTest, InitializeCircularPrimeDimensionCanonicalCounts) {
    // n=5 is prime. For k=1..4, gcd(k,5)=1 so rotational canonicalization applies.
    Supports supports(5, true);
    supports.initialize();

    EXPECT_EQ(supports.get_supports(1).size(), 1u); // C(5,1)/5
    EXPECT_EQ(supports.get_supports(2).size(), 2u); // C(5,2)/5
    EXPECT_EQ(supports.get_supports(3).size(), 2u); // C(5,3)/5
    EXPECT_EQ(supports.get_supports(4).size(), 1u); // C(5,4)/5
    EXPECT_EQ(supports.get_supports(5).size(), 1u); // not coprime, unchanged

    for (size_t k = 1; k <= 4; ++k) {
        for (const auto support : supports.get_supports(k)) {
            EXPECT_TRUE(bs64::is_smallest_representation(support, 5));
        }
    }
}

TEST(SupportsTest, RemoveSupersetsPrunesOnlyLargerBuckets) {
    Supports supports(4, false);
    supports.initialize();

    bitset64 subset = 0ULL;
    subset = bs64::set_bit_at_pos(subset, 0);
    subset = bs64::set_bit_at_pos(subset, 1);

    const size_t before_size_2 = supports.get_supports(2).size();

    supports.remove_supersets(subset, 2);

    // Same-size bucket must not be touched.
    EXPECT_EQ(supports.get_supports(2).size(), before_size_2);

    // Supersets in larger buckets must be removed.
    EXPECT_EQ(supports.get_supports(3).size(), 2u);
    EXPECT_EQ(supports.get_supports(4).size(), 0u);

    for (const auto support : supports.get_supports(3)) {
        EXPECT_FALSE(bs64::is_subset_of(subset, support));
    }
}

TEST(SupportsTest, RemoveSupersetsAutoSizeMatchesExplicitSize) {
    Supports auto_size(5, false);
    Supports explicit_size(5, false);
    auto_size.initialize();
    explicit_size.initialize();

    bitset64 subset = 0ULL;
    subset = bs64::set_bit_at_pos(subset, 0);
    subset = bs64::set_bit_at_pos(subset, 2);
    subset = bs64::set_bit_at_pos(subset, 4);

    auto_size.remove_supersets(subset);
    explicit_size.remove_supersets(subset, bs64::count_set_bits(subset));

    for (size_t k = 1; k <= 5; ++k) {
        EXPECT_EQ(auto_size.get_supports(k), explicit_size.get_supports(k));
    }
}

