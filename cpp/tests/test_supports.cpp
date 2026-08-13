#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/support_generator_circular.hpp>
#include <fracessa/support_generator_non_circular.hpp>

#include "reference_circular_support_generator.hpp"

using namespace fracessa::support;

namespace {

struct StopGeneration {};

bitset_multiword make_multiword_support(size_t dimension, std::initializer_list<size_t> positions)
{
    bitset_multiword result(dimension);
    for (const size_t position : positions) result.set_bit_at_pos(position);
    return result;
}

bitset_multiword make_multiword_support(size_t dimension, bitset support)
{
    bitset_multiword result(dimension);
    while (support != 0) {
        result.set_bit_at_pos(ctz64(support));
        support &= support - 1;
    }
    return result;
}

fracessa::numeric::matrix_int make_circular_integer_matrix(size_t dimension, const std::vector<slong>& half_row)
{
    fracessa::numeric::matrix_int result(dimension, dimension);
    for (size_t row = 0; row < dimension; ++row) {
        for (size_t column = 0; column < dimension; ++column) {
            const size_t forward = (column + dimension - row) % dimension;
            const size_t distance = std::min(forward, dimension - forward);
            if (distance != 0) result(row, column) = fracessa::numeric::integer(half_row[distance - 1]);
        }
    }
    return result;
}

template<class Generator>
std::vector<std::pair<bitset, size_t>> generate_circular(size_t dimension) {
    Generator generator(dimension);
    std::vector<std::pair<bitset, size_t>> generated;
    generator.generate([&](bitset support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
    });
    return generated;
}

template<class Generator>
std::pair<std::vector<std::pair<bitset, size_t>>, size_t> generate_circular_with_forbidden(size_t dimension,
                                                                                           bitset forbidden) {
    Generator generator(dimension);
    std::vector<std::pair<bitset, size_t>> generated;
    size_t multiplier = 0;
    generator.generate([&](bitset support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
        if (support == forbidden)
            multiplier = generator.add_forbidden_orbit(support);
    });
    return {std::move(generated), multiplier};
}

} // namespace

TEST(SupportGeneratorTest, NonCircularMatchesNumericFixedCardinalityOrder) {
    NonCircularSupportGenerator generator(5);
    std::vector<std::pair<bitset, size_t>> generated;
    std::vector<std::pair<bitset, size_t>> expected;

    generator.generate([&](bitset support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
    });

    for (size_t cardinality = 1; cardinality <= 5; ++cardinality) {
        for (bitset support = 1; support < (1ULL << 5); ++support) {
            if (count_set_bits(support) == cardinality)
                expected.emplace_back(support, cardinality);
        }
    }
    EXPECT_EQ(generated, expected);
}

TEST(SupportGeneratorTest, NonCircularPrunesForbiddenSubsetsDuringLaterLayers) {
    NonCircularSupportGenerator generator(5);
    const bitset forbidden = 0b00101;

    std::vector<bitset> generated;
    generator.generate([&](bitset support, size_t cardinality) {
        if (cardinality == 2 && support == forbidden)
            generator.add_forbidden(support);
        if (cardinality == 3)
            generated.push_back(support);
    });

    ASSERT_FALSE(generated.empty());
    for (const bitset support : generated)
        EXPECT_FALSE(is_subset_of(forbidden, support));
}

TEST(SupportGeneratorTest, NonCircularHandlesDimension64) {
    NonCircularSupportGenerator generator(kMaxBitsetDimension);
    size_t generated = 0;
    bitset singleton_union = 0;

    generator.generate([&](bitset support, size_t cardinality) {
        ++generated;
        EXPECT_EQ(cardinality, 1u);
        singleton_union |= support;
        generator.add_forbidden(support);
    });

    EXPECT_EQ(generated, kMaxBitsetDimension);
    EXPECT_EQ(singleton_union, ~bitset{0});
}

TEST(SupportGeneratorTest, MultiwordNonCircularMatchesOneWordOrderAndPruning)
{
    for (size_t dimension = 1; dimension <= 10; ++dimension) {
        NonCircularSupportGenerator one_word(dimension);
        NonCircularSupportGeneratorMultiword multiword(dimension);
        std::vector<std::pair<std::string, size_t>> expected;
        std::vector<std::pair<std::string, size_t>> actual;

        one_word.generate([&](bitset support, size_t cardinality) {
            expected.emplace_back(to_bitstring(support, dimension), cardinality);
            if (dimension == 5 && support == 0b00101) one_word.add_forbidden(support);
        });
        multiword.generate([&](const bitset_multiword& support, size_t cardinality) {
            actual.emplace_back(support.to_bitstring(), cardinality);
            if (dimension == 5 && support == make_multiword_support(dimension, {0, 2})) multiword.add_forbidden(support);
        });

        EXPECT_EQ(actual, expected) << "dimension=" << dimension;
    }
}

TEST(SupportGeneratorTest, MultiwordNonCircularPrunesAfterSingletonsAcrossWordBoundaries)
{
    for (const size_t dimension : {65u, 129u}) {
        NonCircularSupportGeneratorMultiword generator(dimension);
        bitset_multiword singleton_union(dimension);
        size_t generated = 0;

        generator.generate([&](const bitset_multiword& support, size_t cardinality) {
            ++generated;
            EXPECT_EQ(cardinality, 1u);
            singleton_union.union_with(support);
            generator.add_forbidden(support);
        });

        EXPECT_EQ(generated, dimension);
        EXPECT_EQ(singleton_union.count_set_bits(), dimension);
        EXPECT_TRUE(singleton_union.is_set_at_pos(63));
        EXPECT_TRUE(singleton_union.is_set_at_pos(64));
        if (dimension > 128) EXPECT_TRUE(singleton_union.is_set_at_pos(128));
    }
}

TEST(SupportGeneratorTest, MultiwordNonCircularPrunesForbiddenPairAcrossBits63And64)
{
    constexpr size_t dimension = 65;
    NonCircularSupportGeneratorMultiword generator(dimension);
    const bitset_multiword forbidden = make_multiword_support(dimension, {63, 64});
    bool found_forbidden = false;
    size_t surviving_triples = 0;

    EXPECT_THROW(generator.generate([&](const bitset_multiword& support, size_t cardinality) {
        if (cardinality == 2 && support == forbidden) {
            found_forbidden = true;
            generator.add_forbidden(support);
        } else if (cardinality == 3) {
            EXPECT_FALSE(forbidden.is_subset_of(support));
            ++surviving_triples;
        } else if (cardinality == 4) {
            throw StopGeneration{};
        }
    }), StopGeneration);

    EXPECT_TRUE(found_forbidden);
    EXPECT_GT(surviving_triples, 0u);
}

TEST(SupportGeneratorTest, CircularGeneratesEveryNonemptyBraceletOnce) {
    constexpr size_t dimension = 8;
    CircularSupportGenerator generator(dimension);
    using support_table = std::array<bool, 1ULL << dimension>;
    support_table seen{};
    size_t representative_count = 0;
    size_t previous_cardinality = 0;
    bitset previous_representative = 0;

    generator.generate([&](bitset representative, size_t cardinality) {
        ++representative_count;
        EXPECT_EQ(count_set_bits(representative), cardinality);
        if (cardinality == previous_cardinality) {
            EXPECT_GT(representative, previous_representative);
        } else {
            EXPECT_EQ(cardinality, previous_cardinality + 1);
            previous_cardinality = cardinality;
        }
        previous_representative = representative;

        support_table orbit{};
        bitset current = representative;
        for (size_t i = 0; i < dimension; ++i) {
            orbit[current] = true;
            current = rot_one_right(current, dimension);
        }
        current = reflect(representative, dimension);
        for (size_t i = 0; i < dimension; ++i) {
            orbit[current] = true;
            current = rot_one_right(current, dimension);
        }

        for (bitset support = 1; support < orbit.size(); ++support) {
            if (!orbit[support])
                continue;
            EXPECT_FALSE(seen[support]);
            seen[support] = true;
        }
    });

    EXPECT_EQ(representative_count, 29u);
    for (bitset support = 1; support < seen.size(); ++support)
        EXPECT_TRUE(seen[support]);
}

TEST(SupportGeneratorTest, CircularReturnsAndPrunesTheCompleteDihedralOrbit) {
    CircularSupportGenerator generator(6);
    const bitset forbidden = 0b001011;
    const size_t multiplier = generator.add_forbidden_orbit(forbidden);

    EXPECT_EQ(multiplier, 12u);

    std::vector<bitset> generated;
    generator.generate([&](bitset support, size_t cardinality) {
        if (cardinality == 3)
            generated.push_back(support);
    });
    ASSERT_FALSE(generated.empty());
    for (const bitset support : generated) {
        bitset current = forbidden;
        for (size_t i = 0; i < 6; ++i) {
            EXPECT_FALSE(is_subset_of(current, support));
            current = rot_one_right(current, 6);
        }
        current = reflect(forbidden, 6);
        for (size_t i = 0; i < 6; ++i) {
            EXPECT_FALSE(is_subset_of(current, support));
            current = rot_one_right(current, 6);
        }
    }
}

TEST(SupportGeneratorTest, CircularMultiplierCountsDistinctMasksOnly) {
    CircularSupportGenerator singleton_generator(5);
    EXPECT_EQ(singleton_generator.add_forbidden_orbit(0b00001), 5u);

    CircularSupportGenerator full_generator(6);
    EXPECT_EQ(full_generator.add_forbidden_orbit(0b111111), 1u);
}

TEST(SupportGeneratorTest, CircularProductionMatchesReferenceGenerationAndPruning) {
    for (size_t dimension = 1; dimension <= 24; ++dimension)
        EXPECT_EQ(generate_circular<CircularSupportGenerator>(dimension),
                  generate_circular<ReferenceCircularSupportGenerator>(dimension));

    constexpr size_t dimension = 8;
    constexpr bitset forbidden = 0b00000011;
    EXPECT_EQ(generate_circular_with_forbidden<CircularSupportGenerator>(dimension, forbidden),
              generate_circular_with_forbidden<ReferenceCircularSupportGenerator>(dimension, forbidden));
}

TEST(SupportGeneratorTest, CircularHandlesDimension64) {
    constexpr size_t dimension = kMaxBitsetDimension;
    CircularSupportGenerator generator(dimension);
    size_t generated = 0;
    generator.generate([&](bitset support, size_t cardinality) {
        ++generated;
        EXPECT_EQ(cardinality, 1u);
        EXPECT_EQ(support, 1u);
        EXPECT_EQ(generator.add_forbidden_orbit(support), dimension);
    });
    EXPECT_EQ(generated, 1u);
}

TEST(SupportGeneratorTest, MultiwordCircularMatchesOneWordGenerationAndPruning)
{
    for (size_t dimension = 1; dimension <= 16; ++dimension) {
        CircularSupportGenerator one_word(dimension);
        CircularSupportGeneratorMultiword multiword(dimension);
        std::vector<std::pair<std::string, size_t>> expected;
        std::vector<std::pair<std::string, size_t>> actual;

        one_word.generate([&](bitset support, size_t cardinality) {
            expected.emplace_back(to_bitstring(support, dimension), cardinality);
            if (dimension == 8 && support == 0b00000011) one_word.add_forbidden_orbit(support);
        });
        multiword.generate([&](const bitset_multiword& support, size_t cardinality) {
            actual.emplace_back(support.to_bitstring(), cardinality);
            if (dimension == 8 && support == make_multiword_support(dimension, bitset{0b00000011}))
                multiword.add_forbidden_orbit(support);
        });

        EXPECT_EQ(actual, expected) << "dimension=" << dimension;
    }
}

TEST(SupportGeneratorTest, MultiwordCircularPrunesAcrossWordBoundaries)
{
    constexpr size_t dimension = 129;
    CircularSupportGeneratorMultiword generator(dimension);
    const bitset_multiword forbidden = make_multiword_support(dimension, {0, 64});
    std::vector<size_t> indices;
    indices.reserve(dimension);
    bool found_forbidden = false;
    bool crossed_word_boundary = false;
    size_t surviving_triples = 0;

    EXPECT_THROW(generator.generate([&](const bitset_multiword& support, size_t cardinality) {
        EXPECT_EQ(support.count_set_bits(), cardinality);
        if (cardinality == 2 && support == forbidden) {
            found_forbidden = true;
            EXPECT_EQ(generator.add_forbidden_orbit(support), dimension);
        } else if (cardinality == 3) {
            ++surviving_triples;
            support.extract_set_indices(indices);
            crossed_word_boundary = crossed_word_boundary || indices.back() >= 64;

            bitset_multiword forbidden_rotation = forbidden;
            do {
                EXPECT_FALSE(forbidden_rotation.is_subset_of(support));
                forbidden_rotation.rotate_one_right();
            } while (forbidden_rotation != forbidden);
        } else if (cardinality == 4) {
            throw StopGeneration{};
        }
    }), StopGeneration);

    EXPECT_TRUE(found_forbidden);
    EXPECT_TRUE(crossed_word_boundary);
    EXPECT_GT(surviving_triples, 0u);
}

TEST(CircularAffineSymmetryTest, DetectsAndFiltersExactMultiplierSymmetry) {
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(8, {7, 11, 7, 13});
    CircularAffineSymmetry symmetry(matrix);

    ASSERT_EQ(symmetry.multiplier_class_count(), 2u);
    EXPECT_TRUE(symmetry.is_representative(0b00000011));
    EXPECT_FALSE(symmetry.is_representative(0b00001001));

    CircularSupportGenerator generator(8);
    size_t representative_count = 0;
    generator.generate([&](bitset support, size_t) {
        if (symmetry.is_representative(support)) ++representative_count;
    });
    EXPECT_EQ(representative_count, 23u);
}

TEST(CircularAffineSymmetryTest, RepeatedValuesAloneDoNotCreateASymmetry) {
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(8, {7, 7, 11, 13});
    CircularAffineSymmetry symmetry(matrix);

    EXPECT_EQ(symmetry.multiplier_class_count(), 1u);
    EXPECT_TRUE(symmetry.is_representative(0b00001001));
}

TEST(CircularAffineSymmetryTest, EnlargedOrbitReusesExistingDihedralExpansion) {
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(8, {7, 11, 7, 13});
    CircularAffineSymmetry symmetry(matrix);
    CircularSupportGenerator generator(8);

    size_t bracelet_images = 0;
    size_t represented_supports = 0;
    symmetry.for_each_distinct_bracelet_image(0b00000011, [&](bitset image, size_t multiplier_class,
                                                               bool reflected, size_t right_shifts) {
        ++bracelet_images;
        EXPECT_EQ(symmetry.image_mask(0b00000011, multiplier_class, reflected, right_shifts), image);
        const size_t multiplier = generator.add_forbidden_orbit(image);
        EXPECT_LE(multiplier, 16u);
        represented_supports += multiplier;
    });

    EXPECT_EQ(bracelet_images, 2u);
    EXPECT_EQ(represented_supports, 16u);
}

TEST(CircularAffineSymmetryTest, DetectsPublishedDimension24MultiplierClasses) {
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(24, {15, 15, 7, 15, 15, 7, 15, 7, 7, 15, 15, 7});
    CircularAffineSymmetry symmetry(matrix);

    EXPECT_EQ(symmetry.multiplier_class_count(), 4u);
}

TEST(CircularAffineSymmetryTest, MultiwordMatchesOneWordImages)
{
    constexpr size_t dimension = 8;
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(dimension, {7, 11, 7, 13});
    CircularAffineSymmetry one_word(matrix);
    CircularAffineSymmetryMultiword multiword(matrix);

    ASSERT_EQ(multiword.multiplier_class_count(), one_word.multiplier_class_count());
    for (bitset support = 1; support < (bitset{1} << dimension); ++support) {
        const bitset_multiword large_support = make_multiword_support(dimension, support);
        EXPECT_EQ(multiword.is_representative(large_support), one_word.is_representative(support)) << "support=" << support;
    }

    std::vector<std::tuple<std::string, size_t, bool, size_t>> expected;
    std::vector<std::tuple<std::string, size_t, bool, size_t>> actual;
    constexpr bitset support = 0b00000011;
    const bitset_multiword large_support = make_multiword_support(dimension, support);

    one_word.for_each_distinct_bracelet_image(support,
        [&](bitset image, size_t multiplier_class, bool reflected, size_t right_shifts) {
            expected.emplace_back(to_bitstring(image, dimension), multiplier_class, reflected, right_shifts);
        });
    multiword.for_each_distinct_bracelet_image(large_support,
        [&](const bitset_multiword& image, size_t multiplier_class, bool reflected, size_t right_shifts) {
            bitset_multiword reconstructed(dimension);
            multiword.image_mask(large_support, multiplier_class, reflected, right_shifts, reconstructed);
            EXPECT_EQ(reconstructed, image);
            actual.emplace_back(image.to_bitstring(), multiplier_class, reflected, right_shifts);
        });

    EXPECT_EQ(actual, expected);
}

TEST(CircularAffineSymmetryTest, MultiwordImagesCrossWordBoundaries)
{
    for (const size_t dimension : {65u, 129u}) {
        const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(dimension, std::vector<slong>(dimension / 2, 1));
        CircularAffineSymmetryMultiword symmetry(matrix);
        bitset_multiword support(dimension);
        for (const size_t position : {0u, 1u, 63u, 64u, 128u})
            if (position < dimension) support.set_bit_at_pos(position);

        std::vector<size_t> multipliers{1};
        for (size_t multiplier = 2; multiplier <= dimension / 2; ++multiplier)
            if (std::gcd(multiplier, dimension) == 1) multipliers.push_back(multiplier);
        ASSERT_EQ(symmetry.multiplier_class_count(), multipliers.size());

        std::vector<std::string> images;
        symmetry.for_each_distinct_bracelet_image(support,
            [&](const bitset_multiword& image, size_t multiplier_class, bool reflected, size_t right_shifts) {
                bitset_multiword expected(dimension);
                for (const size_t position : {0u, 1u, 63u, 64u, 128u}) {
                    if (position >= dimension) continue;
                    size_t image_position = (multipliers[multiplier_class] * position) % dimension;
                    if (reflected) image_position = dimension - 1 - image_position;
                    image_position = (image_position + dimension - right_shifts) % dimension;
                    expected.set_bit_at_pos(image_position);
                }
                EXPECT_EQ(image, expected);

                bitset_multiword reconstructed(dimension);
                symmetry.image_mask(support, multiplier_class, reflected, right_shifts, reconstructed);
                EXPECT_EQ(reconstructed, image);
                images.push_back(image.to_bitstring());
            });

        EXPECT_GT(images.size(), 1u);
        std::sort(images.begin(), images.end());
        EXPECT_EQ(std::adjacent_find(images.begin(), images.end()), images.end());
    }
}
