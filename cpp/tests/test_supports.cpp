#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <utility>
#include <vector>

#include <fracessa/circular_affine_symmetry.hpp>
#include <fracessa/support_generator_circular.hpp>
#include <fracessa/support_generator_non_circular.hpp>

#include "reference_circular_support_generator.hpp"

using namespace fracessa::support;

namespace {

struct StopGeneration {};

Support make_support(SupportContext& context, bitset bits)
{
    Support result = context.make();
    while (bits != 0) {
        context.set(result, ctz64(bits));
        bits &= bits - 1;
    }
    return result;
}

Support make_support(SupportContext& context, std::initializer_list<size_t> positions)
{
    Support result = context.make();
    for (const size_t position : positions)
        if (position < context.dimension()) context.set(result, position);
    return result;
}

bitset to_bitset(const SupportContext& context, const Support& support)
{
    bitset result = 0;
    for (size_t position = 0; position < context.dimension(); ++position)
        if (context.contains(support, position)) result |= bitset{1} << position;
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

std::vector<std::pair<bitset, size_t>> generate_circular(size_t dimension)
{
    SupportContext context(dimension);
    CircularSupportGenerator generator(context);
    std::vector<std::pair<bitset, size_t>> generated;
    generator.generate([&](const Support& support, size_t cardinality) {
        generated.emplace_back(to_bitset(context, support), cardinality);
    });
    return generated;
}

std::pair<std::vector<std::pair<bitset, size_t>>, size_t> generate_circular_with_forbidden(size_t dimension, bitset forbidden_bits)
{
    SupportContext context(dimension);
    CircularSupportGenerator generator(context);
    const Support forbidden = make_support(context, forbidden_bits);
    std::vector<std::pair<bitset, size_t>> generated;
    size_t multiplier = 0;
    generator.generate([&](const Support& support, size_t cardinality) {
        generated.emplace_back(to_bitset(context, support), cardinality);
        if (context.equal(support, forbidden)) multiplier = generator.add_forbidden_orbit(support);
    });
    return {std::move(generated), multiplier};
}

std::vector<std::pair<bitset, size_t>> generate_reference_circular(size_t dimension)
{
    ReferenceCircularSupportGenerator generator(dimension);
    std::vector<std::pair<bitset, size_t>> generated;
    generator.generate([&](bitset support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
    });
    return generated;
}

std::pair<std::vector<std::pair<bitset, size_t>>, size_t> generate_reference_circular_with_forbidden(size_t dimension,
                                                                                                      bitset forbidden)
{
    ReferenceCircularSupportGenerator generator(dimension);
    std::vector<std::pair<bitset, size_t>> generated;
    size_t multiplier = 0;
    generator.generate([&](bitset support, size_t cardinality) {
        generated.emplace_back(support, cardinality);
        if (support == forbidden) multiplier = generator.add_forbidden_orbit(support);
    });
    return {std::move(generated), multiplier};
}

} // namespace

TEST(SupportGeneratorTest, NonCircularMatchesNumericFixedCardinalityOrder)
{
    SupportContext context(5);
    NonCircularSupportGenerator generator(context);
    std::vector<std::pair<bitset, size_t>> generated;
    std::vector<std::pair<bitset, size_t>> expected;

    generator.generate([&](const Support& support, size_t cardinality) {
        generated.emplace_back(to_bitset(context, support), cardinality);
    });

    for (size_t cardinality = 1; cardinality <= 5; ++cardinality) {
        for (bitset support = 1; support < (1ULL << 5); ++support) {
            if (count_set_bits(support) == cardinality) expected.emplace_back(support, cardinality);
        }
    }
    EXPECT_EQ(generated, expected);
}

TEST(SupportGeneratorTest, NonCircularPrunesForbiddenSubsetsDuringLaterLayers)
{
    SupportContext context(5);
    NonCircularSupportGenerator generator(context);
    const Support forbidden = make_support(context, bitset{0b00101});

    std::vector<bitset> generated;
    generator.generate([&](const Support& support, size_t cardinality) {
        if (cardinality == 2 && context.equal(support, forbidden)) generator.add_forbidden(support);
        if (cardinality == 3) generated.push_back(to_bitset(context, support));
    });

    ASSERT_FALSE(generated.empty());
    for (const bitset support : generated) EXPECT_FALSE(is_subset_of(0b00101, support));
}

TEST(SupportGeneratorTest, NonCircularHandlesDimension64)
{
    SupportContext context(kMaxBitsetDimension);
    NonCircularSupportGenerator generator(context);
    size_t generated = 0;
    bitset singleton_union = 0;

    generator.generate([&](const Support& support, size_t cardinality) {
        ++generated;
        EXPECT_EQ(cardinality, 1u);
        singleton_union |= to_bitset(context, support);
        generator.add_forbidden(support);
    });

    EXPECT_EQ(generated, kMaxBitsetDimension);
    EXPECT_EQ(singleton_union, ~bitset{0});
}

TEST(SupportGeneratorTest, NonCircularReportsWhetherSupportsRemainAfterSingletons)
{
    const auto has_supports_after_forbidding = [](size_t forbidden_count) {
        SupportContext context(4);
        NonCircularSupportGenerator generator(context);
        bool result = false;
        size_t singleton_index = 0;
        EXPECT_THROW(generator.generate([&](const Support& support, size_t cardinality) {
            EXPECT_EQ(cardinality, 1u);
            if (singleton_index++ < forbidden_count) generator.add_forbidden(support);
            if (to_bitset(context, support) == 0b1000) {
                result = generator.has_supports_after_singletons();
                throw StopGeneration{};
            }
        }), StopGeneration);
        return result;
    };

    EXPECT_TRUE(has_supports_after_forbidding(2));
    EXPECT_FALSE(has_supports_after_forbidding(3));
    EXPECT_FALSE(has_supports_after_forbidding(4));
}

TEST(SupportGeneratorTest, NonCircularPrunesAfterSingletonsAcrossWordBoundaries)
{
    for (const size_t dimension : {65u, 129u}) {
        SupportContext context(dimension);
        NonCircularSupportGenerator generator(context);
        Support singleton_union = context.make();
        size_t generated = 0;

        generator.generate([&](const Support& support, size_t cardinality) {
            ++generated;
            EXPECT_EQ(cardinality, 1u);
            context.add(singleton_union, support);
            generator.add_forbidden(support);
            if (generated == dimension) EXPECT_FALSE(generator.has_supports_after_singletons());
        });

        EXPECT_EQ(generated, dimension);
        EXPECT_EQ(context.count(singleton_union), dimension);
        EXPECT_TRUE(context.contains(singleton_union, 63));
        EXPECT_TRUE(context.contains(singleton_union, 64));
        if (dimension > 128) EXPECT_TRUE(context.contains(singleton_union, 128));
    }
}

TEST(SupportGeneratorTest, NonCircularPrunesForbiddenPairAcrossBits63And64)
{
    constexpr size_t dimension = 65;
    SupportContext context(dimension);
    NonCircularSupportGenerator generator(context);
    const Support forbidden = make_support(context, {63, 64});
    bool found_forbidden = false;
    size_t surviving_triples = 0;

    EXPECT_THROW(generator.generate([&](const Support& support, size_t cardinality) {
        if (cardinality == 2 && context.equal(support, forbidden)) {
            found_forbidden = true;
            generator.add_forbidden(support);
        } else if (cardinality == 3) {
            EXPECT_FALSE(context.is_subset_of(forbidden, support));
            ++surviving_triples;
        } else if (cardinality == 4) {
            throw StopGeneration{};
        }
    }), StopGeneration);

    EXPECT_TRUE(found_forbidden);
    EXPECT_GT(surviving_triples, 0u);
}

TEST(SupportGeneratorTest, CircularGeneratesEveryNonemptyBraceletOnce)
{
    constexpr size_t dimension = 8;
    SupportContext context(dimension);
    CircularSupportGenerator generator(context);
    using support_table = std::array<bool, 1ULL << dimension>;
    support_table seen{};
    size_t representative_count = 0;
    size_t previous_cardinality = 0;
    bitset previous_representative = 0;

    generator.generate([&](const Support& support, size_t cardinality) {
        const bitset representative = to_bitset(context, support);
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

        for (bitset orbit_support = 1; orbit_support < orbit.size(); ++orbit_support) {
            if (!orbit[orbit_support]) continue;
            EXPECT_FALSE(seen[orbit_support]);
            seen[orbit_support] = true;
        }
    });

    EXPECT_EQ(representative_count, 29u);
    for (bitset support = 1; support < seen.size(); ++support) EXPECT_TRUE(seen[support]);
}

TEST(SupportGeneratorTest, CircularReturnsAndPrunesTheCompleteDihedralOrbit)
{
    SupportContext context(6);
    CircularSupportGenerator generator(context);
    const Support forbidden = make_support(context, bitset{0b001011});
    const size_t multiplier = generator.add_forbidden_orbit(forbidden);

    EXPECT_EQ(multiplier, 12u);

    std::vector<bitset> generated;
    generator.generate([&](const Support& support, size_t cardinality) {
        if (cardinality == 3) generated.push_back(to_bitset(context, support));
    });
    ASSERT_FALSE(generated.empty());
    for (const bitset support : generated) {
        bitset current = 0b001011;
        for (size_t i = 0; i < 6; ++i) {
            EXPECT_FALSE(is_subset_of(current, support));
            current = rot_one_right(current, 6);
        }
        current = reflect(0b001011, 6);
        for (size_t i = 0; i < 6; ++i) {
            EXPECT_FALSE(is_subset_of(current, support));
            current = rot_one_right(current, 6);
        }
    }
}

TEST(SupportGeneratorTest, CircularMultiplierCountsDistinctMasksOnly)
{
    SupportContext singleton_context(5);
    CircularSupportGenerator singleton_generator(singleton_context);
    const Support singleton = make_support(singleton_context, bitset{0b00001});
    EXPECT_EQ(singleton_generator.add_forbidden_orbit(singleton), 5u);

    SupportContext full_context(6);
    CircularSupportGenerator full_generator(full_context);
    Support full = full_context.make();
    full_context.set_all(full);
    EXPECT_EQ(full_generator.add_forbidden_orbit(full), 1u);
}

TEST(SupportGeneratorTest, CircularProductionMatchesReferenceGenerationAndPruning)
{
    for (size_t dimension = 1; dimension <= 24; ++dimension)
        EXPECT_EQ(generate_circular(dimension), generate_reference_circular(dimension));

    constexpr size_t dimension = 8;
    constexpr bitset forbidden = 0b00000011;
    EXPECT_EQ(generate_circular_with_forbidden(dimension, forbidden),
              generate_reference_circular_with_forbidden(dimension, forbidden));
}

TEST(SupportGeneratorTest, CircularHandlesDimension64)
{
    constexpr size_t dimension = kMaxBitsetDimension;
    SupportContext context(dimension);
    CircularSupportGenerator generator(context);
    size_t generated = 0;
    generator.generate([&](const Support& support, size_t cardinality) {
        ++generated;
        EXPECT_EQ(cardinality, 1u);
        EXPECT_TRUE(context.contains(support, 0));
        EXPECT_TRUE(generator.has_supports_after_singletons());
        EXPECT_EQ(generator.add_forbidden_orbit(support), dimension);
        EXPECT_FALSE(generator.has_supports_after_singletons());
    });
    EXPECT_EQ(generated, 1u);
}

TEST(SupportGeneratorTest, CircularReportsNoSupportAfterDimension65SingletonOrbit)
{
    constexpr size_t dimension = 65;
    SupportContext context(dimension);
    CircularSupportGenerator generator(context);
    size_t generated = 0;

    generator.generate([&](const Support& support, size_t cardinality) {
        ++generated;
        EXPECT_EQ(cardinality, 1u);
        EXPECT_TRUE(generator.has_supports_after_singletons());
        EXPECT_EQ(generator.add_forbidden_orbit(support), dimension);
        EXPECT_FALSE(generator.has_supports_after_singletons());
    });

    EXPECT_EQ(generated, 1u);
}

TEST(SupportGeneratorTest, CircularPrunesAcrossWordBoundaries)
{
    constexpr size_t dimension = 129;
    SupportContext context(dimension);
    CircularSupportGenerator generator(context);
    const Support forbidden = make_support(context, {0, 64});
    std::vector<size_t> indices;
    indices.reserve(dimension);
    bool found_forbidden = false;
    bool crossed_word_boundary = false;
    size_t surviving_triples = 0;

    EXPECT_THROW(generator.generate([&](const Support& support, size_t cardinality) {
        EXPECT_EQ(context.count(support), cardinality);
        if (cardinality == 2 && context.equal(support, forbidden)) {
            found_forbidden = true;
            EXPECT_EQ(generator.add_forbidden_orbit(support), dimension);
        } else if (cardinality == 3) {
            ++surviving_triples;
            context.extract_set_indices(support, indices);
            crossed_word_boundary = crossed_word_boundary || indices.back() >= 64;

            Support forbidden_rotation = context.clone(forbidden);
            do {
                EXPECT_FALSE(context.is_subset_of(forbidden_rotation, support));
                context.rotate_one_right(forbidden_rotation);
            } while (!context.equal(forbidden_rotation, forbidden));
        } else if (cardinality == 4) {
            throw StopGeneration{};
        }
    }), StopGeneration);

    EXPECT_TRUE(found_forbidden);
    EXPECT_TRUE(crossed_word_boundary);
    EXPECT_GT(surviving_triples, 0u);
}

TEST(CircularAffineSymmetryTest, DetectsAndFiltersExactMultiplierSymmetry)
{
    constexpr size_t dimension = 8;
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(dimension, {7, 11, 7, 13});
    SupportContext context(dimension);
    CircularAffineSymmetry symmetry(matrix, context);
    const Support representative = make_support(context, bitset{0b00000011});
    const Support non_representative = make_support(context, bitset{0b00001001});

    ASSERT_EQ(symmetry.multiplier_class_count(), 2u);
    EXPECT_TRUE(symmetry.is_representative(representative));
    EXPECT_FALSE(symmetry.is_representative(non_representative));

    CircularSupportGenerator generator(context);
    size_t representative_count = 0;
    generator.generate([&](const Support& support, size_t) {
        if (symmetry.is_representative(support)) ++representative_count;
    });
    EXPECT_EQ(representative_count, 23u);
}

TEST(CircularAffineSymmetryTest, RepeatedValuesAloneDoNotCreateASymmetry)
{
    constexpr size_t dimension = 8;
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(dimension, {7, 7, 11, 13});
    SupportContext context(dimension);
    CircularAffineSymmetry symmetry(matrix, context);
    const Support support = make_support(context, bitset{0b00001001});

    EXPECT_EQ(symmetry.multiplier_class_count(), 1u);
    EXPECT_TRUE(symmetry.is_representative(support));
}

TEST(CircularAffineSymmetryTest, EnlargedOrbitReusesExistingDihedralExpansion)
{
    constexpr size_t dimension = 8;
    const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(dimension, {7, 11, 7, 13});
    SupportContext context(dimension);
    CircularAffineSymmetry symmetry(matrix, context);
    CircularSupportGenerator generator(context);
    const Support support = make_support(context, bitset{0b00000011});
    Support reconstructed = context.make();

    size_t bracelet_images = 0;
    size_t represented_supports = 0;
    symmetry.for_each_distinct_bracelet_image(support, [&](const Support& image, size_t multiplier_class,
                                                           bool reflected, size_t right_shifts) {
        ++bracelet_images;
        symmetry.image_mask(support, multiplier_class, reflected, right_shifts, reconstructed);
        EXPECT_TRUE(context.equal(reconstructed, image));
        const size_t multiplier = generator.add_forbidden_orbit(image);
        EXPECT_LE(multiplier, 16u);
        represented_supports += multiplier;
    });

    EXPECT_EQ(bracelet_images, 2u);
    EXPECT_EQ(represented_supports, 16u);
}

TEST(CircularAffineSymmetryTest, DetectsPublishedDimension24MultiplierClasses)
{
    constexpr size_t dimension = 24;
    const fracessa::numeric::matrix_int matrix =
        make_circular_integer_matrix(dimension, {15, 15, 7, 15, 15, 7, 15, 7, 7, 15, 15, 7});
    SupportContext context(dimension);
    CircularAffineSymmetry symmetry(matrix, context);

    EXPECT_EQ(symmetry.multiplier_class_count(), 4u);
}

TEST(CircularAffineSymmetryTest, ImagesCrossWordBoundaries)
{
    for (const size_t dimension : {65u, 129u}) {
        const fracessa::numeric::matrix_int matrix = make_circular_integer_matrix(dimension, std::vector<slong>(dimension / 2, 1));
        SupportContext context(dimension);
        CircularAffineSymmetry symmetry(matrix, context);
        const Support support = make_support(context, {0, 1, 63, 64, 128});
        Support expected = context.make();
        Support reconstructed = context.make();

        std::vector<size_t> multipliers{1};
        for (size_t multiplier = 2; multiplier <= dimension / 2; ++multiplier)
            if (std::gcd(multiplier, dimension) == 1) multipliers.push_back(multiplier);
        ASSERT_EQ(symmetry.multiplier_class_count(), multipliers.size());

        std::vector<std::string> images;
        symmetry.for_each_distinct_bracelet_image(support,
            [&](const Support& image, size_t multiplier_class, bool reflected, size_t right_shifts) {
                context.clear(expected);
                for (const size_t position : {0u, 1u, 63u, 64u, 128u}) {
                    if (position >= dimension) continue;
                    size_t image_position = (multipliers[multiplier_class] * position) % dimension;
                    if (reflected) image_position = dimension - 1 - image_position;
                    image_position = (image_position + dimension - right_shifts) % dimension;
                    context.set(expected, image_position);
                }
                EXPECT_TRUE(context.equal(image, expected));

                symmetry.image_mask(support, multiplier_class, reflected, right_shifts, reconstructed);
                EXPECT_TRUE(context.equal(reconstructed, image));
                images.push_back(context.to_bitstring(image));
            });

        EXPECT_GT(images.size(), 1u);
        std::sort(images.begin(), images.end());
        EXPECT_EQ(std::adjacent_find(images.begin(), images.end()), images.end());
    }
}
