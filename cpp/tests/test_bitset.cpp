#include <gtest/gtest.h>
#include <fracessa/bitset.hpp>
#include <vector>

using namespace fracessa::support;

/*
 * Bitset support tests.
 *
 * These assert core set-theoretic primitives used everywhere in support
 * enumeration and pruning. Failures here usually imply broad algorithmic breakage.
 */

// Test Basic Operations
TEST(BitsetTest, Set) {
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 5);
    EXPECT_TRUE(is_set_at_pos(bits, 5));
    EXPECT_FALSE(is_set_at_pos(bits, 4));
}

TEST(BitsetTest, SetAll) {
    bitset bits = 0ULL;
    bits = set_all_n_bits(5);
    EXPECT_EQ(count_set_bits(bits), 5);
    for (unsigned i = 0; i < 5; ++i) {
        EXPECT_TRUE(is_set_at_pos(bits, i));
    }
    EXPECT_FALSE(is_set_at_pos(bits, 5));
}

TEST(BitsetTest, SetAll64) {
    const bitset bits = set_all_n_bits(64);
    EXPECT_EQ(bits, ~bitset{0});
    EXPECT_EQ(count_set_bits(bits), 64u);
}

TEST(BitsetTest, Count) {
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 0);
    bits = set_bit_at_pos(bits, 2);
    bits = set_bit_at_pos(bits, 5);
    EXPECT_EQ(count_set_bits(bits), 3);
}

TEST(BitsetTest, FindFirst) {
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 3);
    EXPECT_EQ(find_pos_first_set_bit(bits), 3);
}

// Test Rotations
TEST(BitsetTest, RotOneRight) {
    bitset bits = 0ULL;
    bits = set_all_n_bits(4); // bits 0,1,2,3 set
    // Original: 1111 (bits 0-3)
    // After rotation: 1111 (bits 0-3) -> 1111 (bits 1-4) -> 1111 (bits 0-3) with wrap

    bits = rot_one_right(bits, 4);
    // After rotation: bit 0 becomes bit 3, bit 1 becomes bit 0, etc.
    // So 1111 -> 1111 (rotated)
    EXPECT_EQ(count_set_bits(bits), 4);
}

TEST(BitsetTest, RotOneRightExactPattern) {
    // n=5, start with bits {0,2} -> bitstring "00101"
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 0);
    bits = set_bit_at_pos(bits, 2);

    bits = rot_one_right(bits, 5);

    // After right-rotate by one on [0..4], bits become {1,4} -> "10010".
    EXPECT_TRUE(is_set_at_pos(bits, 1));
    EXPECT_TRUE(is_set_at_pos(bits, 4));
    EXPECT_FALSE(is_set_at_pos(bits, 0));
    EXPECT_FALSE(is_set_at_pos(bits, 2));
    EXPECT_EQ(to_bitstring(bits, 5), "10010");
}

TEST(BitsetTest, RotOneRightMasksBitsOutsideDimension) {
    // Bit 6 is outside n=5 and must be masked out.
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 0);
    bits = set_bit_at_pos(bits, 6);

    bits = rot_one_right(bits, 5);

    EXPECT_EQ(count_set_bits(bits), 1u);
    EXPECT_TRUE(is_set_at_pos(bits, 4));
    EXPECT_FALSE(is_set_at_pos(bits, 6));
}

TEST(BitsetTest, ReflectExactPatternAndRoundTrip) {
    const bitset bits = 0b001011;
    const bitset reflected = reflect(bits, 6);

    EXPECT_EQ(reflected, 0b110100u);
    EXPECT_EQ(reflect(reflected, 6), bits);
}

TEST(BitsetTest, ReflectMatchesDirectDefinitionForEveryDimension) {
    constexpr bitset samples[] = {
        0ULL,
        ~0ULL,
        0x0123456789abcdefULL,
        0xaaaaaaaaaaaaaaaaULL,
        1ULL << 63,
    };

    for (size_t dimension = 0; dimension <= 64; ++dimension) {
        for (const bitset bits : samples) {
            bitset expected = 0;
            for (size_t i = 0; i < dimension; ++i)
                expected |= ((bits >> i) & 1ULL) << (dimension - 1 - i);
            EXPECT_EQ(reflect(bits, dimension), expected) << "dimension=" << dimension;
        }
    }
}

TEST(BitsetTest, RotationsHandleDimension64) {
    EXPECT_EQ(rot_one_right(1ULL, 64), 1ULL << 63);
}


// Test Subset Operations
TEST(BitsetTest, IsSubsetOf) {
    bitset a = 0ULL;
    a = set_bit_at_pos(a, 1);
    a = set_bit_at_pos(a, 2);

    bitset b = 0ULL;
    b = set_bit_at_pos(b, 0);
    b = set_bit_at_pos(b, 1);
    b = set_bit_at_pos(b, 2);
    b = set_bit_at_pos(b, 3);

    EXPECT_TRUE(is_subset_of(a, b));
    EXPECT_FALSE(is_subset_of(b, a));
}

TEST(BitsetTest, Subtract) {
    bitset a = 0ULL;
    a = set_bit_at_pos(a, 1);
    a = set_bit_at_pos(a, 2);
    a = set_bit_at_pos(a, 3);

    bitset b = 0ULL;
    b = set_bit_at_pos(b, 2);

    bitset result = subtract(a, b);
    EXPECT_TRUE(is_set_at_pos(result, 1));
    EXPECT_FALSE(is_set_at_pos(result, 2));
    EXPECT_TRUE(is_set_at_pos(result, 3));
}

TEST(BitsetTest, LowestSetBit) {
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 5);
    bits = set_bit_at_pos(bits, 7);

    bitset lowest = lowest_set_bit_as_bit(bits);
    EXPECT_EQ(find_pos_first_set_bit(lowest), 5);
    EXPECT_EQ(count_set_bits(lowest), 1);
}

TEST(BitsetTest, LowestSetBitZero) {
    bitset bits = 0ULL;
    bitset lowest = lowest_set_bit_as_bit(bits);
    EXPECT_EQ(lowest, 0ULL);
}

TEST(BitsetTest, GosperNextSameCardinality) {
    EXPECT_EQ(next_same_cardinality(0b00111), 0b01011u);
    EXPECT_EQ(next_same_cardinality(0b01110), 0b10011u);
}

TEST(BitsetTest, ExtractSetIndicesBasic) {
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 0);
    bits = set_bit_at_pos(bits, 2);
    bits = set_bit_at_pos(bits, 5);

    uint8_t idx[kMaxBitsetDimension] = {};
    const size_t count = extract_set_indices(bits, 6, idx);

    ASSERT_EQ(count, 3u);
    EXPECT_EQ(idx[0], 0u);
    EXPECT_EQ(idx[1], 2u);
    EXPECT_EQ(idx[2], 5u);
}

TEST(BitsetTest, ExtractSetIndicesRespectsDimension) {
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 1);
    bits = set_bit_at_pos(bits, 3);
    bits = set_bit_at_pos(bits, 6);

    uint8_t idx[kMaxBitsetDimension] = {};
    const size_t count = extract_set_indices(bits, 5, idx);

    ASSERT_EQ(count, 2u);
    EXPECT_EQ(idx[0], 1u);
    EXPECT_EQ(idx[1], 3u);
}

TEST(BitsetTest, ExtractSetIndicesHandlesTopBit) {
    uint8_t indices[kMaxBitsetDimension] = {};
    const size_t count = extract_set_indices(1ULL | (1ULL << 63), 64, indices);

    ASSERT_EQ(count, 2u);
    EXPECT_EQ(indices[0], 0u);
    EXPECT_EQ(indices[1], 63u);
}

// Test String Functions
TEST(BitsetTest, ToString) {
    bitset bits = 42ULL;
    std::string s = to_string(bits);
    EXPECT_EQ(s, "42");
}

TEST(BitsetTest, ToBitstring) {
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 0);
    bits = set_bit_at_pos(bits, 2);
    bits = set_bit_at_pos(bits, 3);
    // bits 0,2,3 set in dimension 4
    // MSB first: bit 3, bit 2, bit 1, bit 0
    // So: 1101

    std::string s = to_bitstring(bits, 4);
    EXPECT_EQ(s, "1101");
}

// Test Iterate All Supports
TEST(BitsetTest, IterateAllSupports) {
    std::vector<bitset> supports;
    for (bitset bits = 1ull; bits < (1ull << 3); ++bits) {
        supports.push_back(bits);
    }
    // Should generate 2^3 - 1 = 7 non-empty supports
    EXPECT_EQ(supports.size(), 7);
}

// Test Edge Cases
TEST(BitsetTest, EmptyBitset) {
    bitset bits = 0ULL;
    EXPECT_EQ(count_set_bits(bits), 0);
}

TEST(BitsetTest, FullBitset) {
    bitset bits = 0ULL;
    bits = set_all_n_bits(8);
    EXPECT_EQ(count_set_bits(bits), 8);
    for (unsigned i = 0; i < 8; ++i) {
        EXPECT_TRUE(is_set_at_pos(bits, i));
    }
}

TEST(BitsetTest, SetAllZero) {
    bitset bits = 0ULL;
    bits = set_all_n_bits(0);
    EXPECT_EQ(count_set_bits(bits), 0);
}
