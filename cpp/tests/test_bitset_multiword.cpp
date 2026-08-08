#include <gtest/gtest.h>

#include <fracessa/bitset_multiword.hpp>

#include <algorithm>
#include <initializer_list>
#include <string>
#include <vector>

namespace {

std::vector<bool> make_reference(size_t dimension, std::initializer_list<size_t> positions)
{
    std::vector<bool> result(dimension, false);
    for (const size_t position : positions)
        if (position < dimension) result[position] = true;
    return result;
}

bitset_multiword make_bitset(const std::vector<bool>& reference)
{
    bitset_multiword result(reference.size());
    for (size_t position = 0; position < reference.size(); ++position)
        if (reference[position]) result.set_bit_at_pos(position);
    return result;
}

std::vector<size_t> reference_indices(const std::vector<bool>& reference)
{
    std::vector<size_t> result;
    for (size_t position = 0; position < reference.size(); ++position)
        if (reference[position]) result.push_back(position);
    return result;
}

std::string reference_bitstring(const std::vector<bool>& reference)
{
    std::string result;
    result.reserve(reference.size());
    for (size_t position = reference.size(); position > 0; --position)
        result += reference[position - 1] ? '1' : '0';
    return result;
}

void expect_matches(const bitset_multiword& actual, const std::vector<bool>& reference)
{
    EXPECT_EQ(actual.dimension(), reference.size());
    EXPECT_EQ(actual.empty(), std::none_of(reference.begin(), reference.end(), [](bool value) { return value; }));
    EXPECT_EQ(actual.count_set_bits(), static_cast<size_t>(std::count(reference.begin(), reference.end(), true)));
    EXPECT_EQ(actual.to_bitstring(), reference_bitstring(reference));

    std::vector<size_t> indices;
    indices.reserve(reference.size());
    const size_t capacity = indices.capacity();
    actual.extract_set_indices(indices);
    EXPECT_EQ(indices, reference_indices(reference));
    EXPECT_EQ(indices.capacity(), capacity);

    std::vector<size_t> unset_indices;
    unset_indices.reserve(reference.size());
    const size_t unset_capacity = unset_indices.capacity();
    actual.extract_unset_indices(unset_indices);
    std::vector<size_t> expected_unset_indices;
    for (size_t position = 0; position < reference.size(); ++position)
        if (!reference[position]) expected_unset_indices.push_back(position);
    EXPECT_EQ(unset_indices, expected_unset_indices);
    EXPECT_EQ(unset_indices.capacity(), unset_capacity);

    for (size_t position = 0; position < reference.size(); ++position)
        EXPECT_EQ(actual.is_set_at_pos(position), reference[position]);
}

void rotate_reference_one_right(std::vector<bool>& reference)
{
    const bool wrap = reference[0];
    for (size_t position = 0; position + 1 < reference.size(); ++position) reference[position] = reference[position + 1];
    reference.back() = wrap;
}

} // namespace

TEST(BitsetMultiwordTest, RejectsZeroDimension)
{
    EXPECT_THROW(bitset_multiword(0), std::invalid_argument);
}

TEST(BitsetMultiwordTest, SetClearAndFillMatchReferenceAcrossWordBoundaries)
{
    for (const size_t dimension : {65u, 70u, 128u, 129u}) {
        bitset_multiword bits(dimension);
        std::vector<bool> reference(dimension, false);
        EXPECT_EQ(bits.word_count(), dimension / 64 + static_cast<size_t>(dimension % 64 != 0));
        expect_matches(bits, reference);

        for (const size_t position : {0u, 63u, 64u, 127u, 128u}) {
            if (position >= dimension) continue;
            bits.set_bit_at_pos(position);
            reference[position] = true;
        }
        expect_matches(bits, reference);

        bits.clear_bit_at_pos(63);
        reference[63] = false;
        expect_matches(bits, reference);

        bits.set_all();
        std::fill(reference.begin(), reference.end(), true);
        expect_matches(bits, reference);

        bits.clear_all();
        std::fill(reference.begin(), reference.end(), false);
        expect_matches(bits, reference);
    }
}

TEST(BitsetMultiwordTest, SetAlgebraMatchesReference)
{
    const std::vector<bool> left_reference = make_reference(129, {0, 63, 64, 128});
    const std::vector<bool> right_reference = make_reference(129, {0, 64, 100, 128});
    bitset_multiword left = make_bitset(left_reference);
    const bitset_multiword right = make_bitset(right_reference);

    bitset_multiword difference = left;
    difference.subtract(right);
    expect_matches(difference, make_reference(129, {63}));

    left.union_with(right);
    expect_matches(left, make_reference(129, {0, 63, 64, 100, 128}));
    EXPECT_TRUE(right.is_subset_of(left));
    EXPECT_FALSE(left.is_subset_of(right));
}

TEST(BitsetMultiwordTest, CopiesIntoExistingStorage)
{
    const bitset_multiword source = make_bitset(make_reference(129, {0, 63, 64, 128}));
    bitset_multiword destination(129);

    destination.copy_from(source);

    EXPECT_EQ(destination, source);
    EXPECT_EQ(destination.word_count(), 3u);
}

TEST(BitsetMultiwordTest, OrderingUsesTheHighestDifferentWord)
{
    const bitset_multiword bit_63 = make_bitset(make_reference(129, {63}));
    const bitset_multiword bit_64 = make_bitset(make_reference(129, {64}));
    const bitset_multiword bit_128 = make_bitset(make_reference(129, {128}));

    EXPECT_LT(bit_63, bit_64);
    EXPECT_LT(bit_64, bit_128);
    EXPECT_FALSE(bit_128 < bit_128);
    EXPECT_EQ(bit_64, make_bitset(make_reference(129, {64})));
}

TEST(BitsetMultiwordTest, FirstBitTraversalAndExtractionCrossWordBoundaries)
{
    bitset_multiword bits = make_bitset(make_reference(129, {0, 63, 64, 127, 128}));
    const std::vector<size_t> expected{0, 63, 64, 127, 128};

    for (const size_t position : expected) {
        ASSERT_FALSE(bits.empty());
        EXPECT_EQ(bits.find_pos_first_set_bit(), position);
        bits.clear_bit_at_pos(position);
    }
    EXPECT_TRUE(bits.empty());
}

TEST(BitsetMultiwordTest, RotationMatchesReferenceAndReturnsAfterOneCycle)
{
    for (const size_t dimension : {65u, 70u, 128u, 129u}) {
        std::vector<bool> reference = make_reference(dimension, {0, 1, 63, 64, 69, 127, 128});
        bitset_multiword bits = make_bitset(reference);
        const bitset_multiword original = bits;

        for (size_t rotation = 0; rotation < dimension; ++rotation) {
            bits.rotate_one_right();
            rotate_reference_one_right(reference);
            expect_matches(bits, reference);
        }
        EXPECT_EQ(bits, original);
    }
}

TEST(BitsetMultiwordTest, ReflectionMatchesReferenceAndIsAnInvolution)
{
    for (size_t dimension = 1; dimension <= 260; ++dimension) {
        std::vector<bool> original_reference(dimension, false);
        for (size_t position = 0; position < dimension; ++position)
            original_reference[position] = ((position * 37 + dimension * 13) % 11) < 5;
        std::vector<bool> reflected_reference(original_reference.rbegin(), original_reference.rend());
        bitset_multiword bits = make_bitset(original_reference);

        bits.reflect();
        expect_matches(bits, reflected_reference);

        bits.reflect();
        expect_matches(bits, original_reference);
    }
}
