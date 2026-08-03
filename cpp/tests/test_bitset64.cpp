#include <gtest/gtest.h>
#include <fracessa/bitset64.hpp>
#include <vector>

/*
 * Bitset support tests.
 *
 * These assert core set-theoretic primitives used everywhere in support
 * enumeration and pruning. Failures here usually imply broad algorithmic breakage.
 */

// Test Basic Operations
TEST(Bitset64Test, Set) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 5);
    EXPECT_TRUE(bs64::is_set_at_pos(bits, 5));
    EXPECT_FALSE(bs64::is_set_at_pos(bits, 4));
}

TEST(Bitset64Test, SetAll) {
    bitset64 bits = 0ULL;
    bits = bs64::set_all_n_bits(5);
    EXPECT_EQ(bs64::count_set_bits(bits), 5);
    for (unsigned i = 0; i < 5; ++i) {
        EXPECT_TRUE(bs64::is_set_at_pos(bits, i));
    }
    EXPECT_FALSE(bs64::is_set_at_pos(bits, 5));
}

TEST(Bitset64Test, Count) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 2);
    bits = bs64::set_bit_at_pos(bits, 5);
    EXPECT_EQ(bs64::count_set_bits(bits), 3);
}

TEST(Bitset64Test, FindFirst) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 3);
    EXPECT_EQ(bs64::find_pos_first_set_bit(bits), 3);
}

TEST(Bitset64Test, FindNext) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 2);
    bits = bs64::set_bit_at_pos(bits, 5);
    bits = bs64::set_bit_at_pos(bits, 7);
    
    unsigned pos = bs64::find_pos_first_set_bit(bits);
    EXPECT_EQ(pos, 2);
    
    pos = bs64::find_pos_next_set_bit(bits, pos);
    EXPECT_EQ(pos, 5);
    
    pos = bs64::find_pos_next_set_bit(bits, pos);
    EXPECT_EQ(pos, 7);
    
    pos = bs64::find_pos_next_set_bit(bits, pos);
    EXPECT_EQ(pos, 64); // No more bits
}

// Test Iteration
TEST(Bitset64Test, ForEachSetBit) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 1);
    bits = bs64::set_bit_at_pos(bits, 3);
    bits = bs64::set_bit_at_pos(bits, 5);
    
    std::vector<unsigned> positions;
    for (unsigned i = bs64::find_pos_first_set_bit(bits); i < 64; i = bs64::find_pos_next_set_bit(bits, i)) {
        positions.push_back(i);
    }
    
    EXPECT_EQ(positions.size(), 3);
    EXPECT_EQ(positions[0], 1);
    EXPECT_EQ(positions[1], 3);
    EXPECT_EQ(positions[2], 5);
}

TEST(Bitset64Test, ForEachSetBitEarlyExit) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 1);
    bits = bs64::set_bit_at_pos(bits, 3);
    bits = bs64::set_bit_at_pos(bits, 5);
    
    std::vector<unsigned> positions;
    for (unsigned i = bs64::find_pos_first_set_bit(bits); i < 64; i = bs64::find_pos_next_set_bit(bits, i)) {
        positions.push_back(i);
        break; // Early exit
    }
    
    EXPECT_EQ(positions.size(), 1);
}

TEST(Bitset64Test, ForEachSetBitNoExit) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 1);
    bits = bs64::set_bit_at_pos(bits, 3);
    bits = bs64::set_bit_at_pos(bits, 5);
    
    std::vector<unsigned> positions;
    for (unsigned i = bs64::find_pos_first_set_bit(bits); i < 64; i = bs64::find_pos_next_set_bit(bits, i)) {
        positions.push_back(i);
    }
    
    EXPECT_EQ(positions.size(), 3);
}

// Test Rotations
TEST(Bitset64Test, RotOneRight) {
    bitset64 bits = 0ULL;
    bits = bs64::set_all_n_bits(4); // bits 0,1,2,3 set
    // Original: 1111 (bits 0-3)
    // After rotation: 1111 (bits 0-3) -> 1111 (bits 1-4) -> 1111 (bits 0-3) with wrap
    
    bits = bs64::rot_one_right(bits, 4);
    // After rotation: bit 0 becomes bit 3, bit 1 becomes bit 0, etc.
    // So 1111 -> 1111 (rotated)
    EXPECT_EQ(bs64::count_set_bits(bits), 4);
}

TEST(Bitset64Test, RotOneRightExactPattern) {
    // n=5, start with bits {0,2} -> bitstring "00101"
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 2);

    bits = bs64::rot_one_right(bits, 5);

    // After right-rotate by one on [0..4], bits become {1,4} -> "10010".
    EXPECT_TRUE(bs64::is_set_at_pos(bits, 1));
    EXPECT_TRUE(bs64::is_set_at_pos(bits, 4));
    EXPECT_FALSE(bs64::is_set_at_pos(bits, 0));
    EXPECT_FALSE(bs64::is_set_at_pos(bits, 2));
    EXPECT_EQ(bs64::to_bitstring(bits, 5), "10010");
}

TEST(Bitset64Test, RotOneRightMasksBitsOutsideDimension) {
    // Bit 6 is outside n=5 and must be masked out.
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 6);

    bits = bs64::rot_one_right(bits, 5);

    EXPECT_EQ(bs64::count_set_bits(bits), 1u);
    EXPECT_TRUE(bs64::is_set_at_pos(bits, 4));
    EXPECT_FALSE(bs64::is_set_at_pos(bits, 6));
}

TEST(Bitset64Test, ReflectExactPatternAndRoundTrip) {
    const bitset64 bits = 0b001011;
    const bitset64 reflected = bs64::reflect(bits, 6);

    EXPECT_EQ(reflected, 0b110100u);
    EXPECT_EQ(bs64::reflect(reflected, 6), bits);
}

TEST(Bitset64Test, ReflectMatchesDirectDefinitionForEveryDimension) {
    constexpr bitset64 samples[] = {
        0ULL,
        ~0ULL,
        0x0123456789abcdefULL,
        0xaaaaaaaaaaaaaaaaULL,
        1ULL << 63,
    };

    for (size_t dimension = 0; dimension <= 64; ++dimension) {
        for (const bitset64 bits : samples) {
            bitset64 expected = 0;
            for (size_t i = 0; i < dimension; ++i)
                expected |= ((bits >> i) & 1ULL) << (dimension - 1 - i);
            EXPECT_EQ(bs64::reflect(bits, dimension), expected) << "dimension=" << dimension;
        }
    }
}


TEST(Bitset64Test, IsSmallestRepresentation) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 1);
    // bits = 0011 (bits 0-3) - this is the smallest representation
    
    EXPECT_TRUE(bs64::is_smallest_representation(bits, 4));
}

// Test Subset Operations
TEST(Bitset64Test, IsSubsetOf) {
    bitset64 a = 0ULL;
    a = bs64::set_bit_at_pos(a, 1);
    a = bs64::set_bit_at_pos(a, 2);
    
    bitset64 b = 0ULL;
    b = bs64::set_bit_at_pos(b, 0);
    b = bs64::set_bit_at_pos(b, 1);
    b = bs64::set_bit_at_pos(b, 2);
    b = bs64::set_bit_at_pos(b, 3);
    
    EXPECT_TRUE(bs64::is_subset_of(a, b));
    EXPECT_FALSE(bs64::is_subset_of(b, a));
}

TEST(Bitset64Test, Subtract) {
    bitset64 a = 0ULL;
    a = bs64::set_bit_at_pos(a, 1);
    a = bs64::set_bit_at_pos(a, 2);
    a = bs64::set_bit_at_pos(a, 3);
    
    bitset64 b = 0ULL;
    b = bs64::set_bit_at_pos(b, 2);
    
    bitset64 result = bs64::subtract(a, b);
    EXPECT_TRUE(bs64::is_set_at_pos(result, 1));
    EXPECT_FALSE(bs64::is_set_at_pos(result, 2));
    EXPECT_TRUE(bs64::is_set_at_pos(result, 3));
}

TEST(Bitset64Test, LowestSetBit) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 5);
    bits = bs64::set_bit_at_pos(bits, 7);
    
    bitset64 lowest = bs64::lowest_set_bit_as_bit(bits);
    EXPECT_EQ(bs64::find_pos_first_set_bit(lowest), 5);
    EXPECT_EQ(bs64::count_set_bits(lowest), 1);
}

TEST(Bitset64Test, LowestSetBitZero) {
    bitset64 bits = 0ULL;
    bitset64 lowest = bs64::lowest_set_bit_as_bit(bits);
    EXPECT_EQ(lowest, 0ULL);
}

TEST(Bitset64Test, GosperNextSameCardinality) {
    EXPECT_EQ(bs64::next_same_cardinality(0b00111), 0b01011u);
    EXPECT_EQ(bs64::next_same_cardinality(0b01110), 0b10011u);
}

TEST(Bitset64Test, FindNextSetBitAtTopBit) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 63);

    const size_t first = bs64::find_pos_first_set_bit(bits);
    EXPECT_EQ(first, 63u);
    EXPECT_EQ(bs64::find_pos_next_set_bit(bits, first), 64u);
}

TEST(Bitset64Test, BitsBeforePosBoundary) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 3);
    bits = bs64::set_bit_at_pos(bits, 5);

    EXPECT_EQ(bs64::bits_before_pos(bits, 0), 0ULL);
    EXPECT_EQ(bs64::bits_before_pos(bits, 4), (1ULL << 0) | (1ULL << 3));
}

TEST(Bitset64Test, ExtractSetIndicesBasic) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 2);
    bits = bs64::set_bit_at_pos(bits, 5);

    uint8_t idx[bs64::kMaxBitsetDimension] = {};
    const size_t count = bs64::extract_set_indices(bits, 6, idx);

    ASSERT_EQ(count, 3u);
    EXPECT_EQ(idx[0], 0u);
    EXPECT_EQ(idx[1], 2u);
    EXPECT_EQ(idx[2], 5u);
}

TEST(Bitset64Test, ExtractSetIndicesRespectsDimension) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 1);
    bits = bs64::set_bit_at_pos(bits, 3);
    bits = bs64::set_bit_at_pos(bits, 6);

    uint8_t idx[bs64::kMaxBitsetDimension] = {};
    const size_t count = bs64::extract_set_indices(bits, 5, idx);

    ASSERT_EQ(count, 2u);
    EXPECT_EQ(idx[0], 1u);
    EXPECT_EQ(idx[1], 3u);
}

// Test String Functions
TEST(Bitset64Test, ToString) {
    bitset64 bits = 42ULL;
    std::string s = bs64::to_string(bits);
    EXPECT_EQ(s, "42");
}

TEST(Bitset64Test, ToBitstring) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 2);
    bits = bs64::set_bit_at_pos(bits, 3);
    // bits 0,2,3 set in dimension 4
    // MSB first: bit 3, bit 2, bit 1, bit 0
    // So: 1101
    
    std::string s = bs64::to_bitstring(bits, 4);
    EXPECT_EQ(s, "1101");
}

// Test Iterate All Supports
TEST(Bitset64Test, IterateAllSupports) {
    std::vector<bitset64> supports;
    for (bitset64 bits = 1ull; bits < (1ull << 3); ++bits) {
        supports.push_back(bits);
    }
    // Should generate 2^3 - 1 = 7 non-empty supports
    EXPECT_EQ(supports.size(), 7);
}

// Test Edge Cases
TEST(Bitset64Test, EmptyBitset) {
    bitset64 bits = 0ULL;
    EXPECT_EQ(bs64::count_set_bits(bits), 0);
}

TEST(Bitset64Test, FullBitset) {
    bitset64 bits = 0ULL;
    bits = bs64::set_all_n_bits(8);
    EXPECT_EQ(bs64::count_set_bits(bits), 8);
    for (unsigned i = 0; i < 8; ++i) {
        EXPECT_TRUE(bs64::is_set_at_pos(bits, i));
    }
}

TEST(Bitset64Test, SetAllZero) {
    bitset64 bits = 0ULL;
    bits = bs64::set_all_n_bits(0);
    EXPECT_EQ(bs64::count_set_bits(bits), 0);
}
