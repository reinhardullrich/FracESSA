#include <gtest/gtest.h>

#include <array>
#include <utility>
#include <vector>

#include <fracessa/supports.hpp>

namespace {

template<class Generator>
std::vector<std::pair<bitset64, size_t>> generate_circular(size_t dimension) {
    Generator generator(dimension);
    std::vector<std::pair<bitset64, size_t>> generated;
    generator.generate([&](bitset64 support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
    });
    return generated;
}

template<class Generator>
std::pair<std::vector<std::pair<bitset64, size_t>>, size_t> generate_circular_with_forbidden(size_t dimension,
                                                                                           bitset64 forbidden) {
    Generator generator(dimension);
    std::vector<std::pair<bitset64, size_t>> generated;
    size_t multiplier = 0;
    generator.generate([&](bitset64 support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
        if (support == forbidden)
            multiplier = generator.add_forbidden(support);
    });
    return {std::move(generated), multiplier};
}

} // namespace

TEST(SupportGeneratorTest, NonCircularMatchesNumericFixedCardinalityOrder) {
    NonCircularSupportGenerator generator(5);
    std::vector<std::pair<bitset64, size_t>> generated;
    std::vector<std::pair<bitset64, size_t>> expected;

    generator.generate([&](bitset64 support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
    });

    for (size_t cardinality = 1; cardinality <= 5; ++cardinality) {
        for (bitset64 support = 1; support < (1ULL << 5); ++support) {
            if (bs64::count_set_bits(support) == cardinality)
                expected.emplace_back(support, cardinality);
        }
    }
    EXPECT_EQ(generated, expected);
}

TEST(SupportGeneratorTest, NonCircularPrunesForbiddenSubsetsDuringLaterLayers) {
    NonCircularSupportGenerator generator(5);
    const bitset64 forbidden = 0b00101;

    std::vector<bitset64> generated;
    generator.generate([&](bitset64 support, size_t cardinality) {
        if (cardinality == 2 && support == forbidden)
            generator.add_forbidden(support);
        if (cardinality == 3)
            generated.push_back(support);
    });

    ASSERT_FALSE(generated.empty());
    for (const bitset64 support : generated)
        EXPECT_FALSE(bs64::is_subset_of(forbidden, support));
}

TEST(SupportGeneratorTest, CircularGeneratesEveryNonemptyBraceletOnce) {
    constexpr size_t dimension = 8;
    CircularSupportGenerator generator(dimension);
    std::array<bool, 1ULL << dimension> seen{};
    size_t representative_count = 0;
    size_t previous_cardinality = 0;
    bitset64 previous_representative = 0;

    generator.generate([&](bitset64 representative, size_t cardinality) {
        ++representative_count;
        EXPECT_EQ(bs64::count_set_bits(representative), cardinality);
        if (cardinality == previous_cardinality) {
            EXPECT_GT(representative, previous_representative);
        } else {
            EXPECT_EQ(cardinality, previous_cardinality + 1);
            previous_cardinality = cardinality;
        }
        previous_representative = representative;

        std::array<bool, 1ULL << dimension> orbit{};
        bitset64 current = representative;
        for (size_t i = 0; i < dimension; ++i) {
            orbit[current] = true;
            current = bs64::rot_one_right(current, dimension);
        }
        current = bs64::reflect(representative, dimension);
        for (size_t i = 0; i < dimension; ++i) {
            orbit[current] = true;
            current = bs64::rot_one_right(current, dimension);
        }

        for (bitset64 support = 1; support < orbit.size(); ++support) {
            if (!orbit[support])
                continue;
            EXPECT_FALSE(seen[support]);
            seen[support] = true;
        }
    });

    EXPECT_EQ(representative_count, 29u);
    for (bitset64 support = 1; support < seen.size(); ++support)
        EXPECT_TRUE(seen[support]);
}

TEST(SupportGeneratorTest, CircularReturnsAndPrunesTheCompleteDihedralOrbit) {
    CircularSupportGenerator generator(6);
    const bitset64 forbidden = 0b001011;
    const size_t multiplier = generator.add_forbidden(forbidden);

    EXPECT_EQ(multiplier, 12u);

    std::vector<bitset64> generated;
    generator.generate([&](bitset64 support, size_t cardinality) {
        if (cardinality == 3)
            generated.push_back(support);
    });
    ASSERT_FALSE(generated.empty());
    for (const bitset64 support : generated) {
        bitset64 current = forbidden;
        for (size_t i = 0; i < 6; ++i) {
            EXPECT_FALSE(bs64::is_subset_of(current, support));
            current = bs64::rot_one_right(current, 6);
        }
        current = bs64::reflect(forbidden, 6);
        for (size_t i = 0; i < 6; ++i) {
            EXPECT_FALSE(bs64::is_subset_of(current, support));
            current = bs64::rot_one_right(current, 6);
        }
    }
}

TEST(SupportGeneratorTest, CircularMultiplierCountsDistinctMasksOnly) {
    CircularSupportGenerator singleton_generator(5);
    EXPECT_EQ(singleton_generator.add_forbidden(0b00001), 5u);

    CircularSupportGenerator full_generator(6);
    EXPECT_EQ(full_generator.add_forbidden(0b111111), 1u);
}

TEST(SupportGeneratorTest, CircularV3MatchesV1GenerationAndPruning) {
    for (size_t dimension = 1; dimension <= 24; ++dimension)
        EXPECT_EQ(generate_circular<CircularSupportGeneratorV3>(dimension),
                  generate_circular<CircularSupportGenerator>(dimension));

    constexpr size_t dimension = 8;
    constexpr bitset64 forbidden = 0b00000011;
    EXPECT_EQ(generate_circular_with_forbidden<CircularSupportGeneratorV3>(dimension, forbidden),
              generate_circular_with_forbidden<CircularSupportGenerator>(dimension, forbidden));
}

TEST(SupportGeneratorTest, CircularV3HandlesDimension63) {
    constexpr size_t dimension = bs64::kMaxBitsetDimension - 1;
    CircularSupportGeneratorV3 generator(dimension);
    size_t generated = 0;
    generator.generate([&](bitset64 support, size_t cardinality) {
        ++generated;
        EXPECT_EQ(cardinality, 1u);
        EXPECT_EQ(support, 1u);
        EXPECT_EQ(generator.add_forbidden(support), dimension);
    });
    EXPECT_EQ(generated, 1u);
}
