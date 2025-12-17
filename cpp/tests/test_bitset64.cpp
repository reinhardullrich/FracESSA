#include <gtest/gtest.h>
#include <fracessa/bitset64.hpp>

// Test Basic Operations
TEST(Bitset64Test, Set) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 5);
    EXPECT_TRUE(bs64::is_set_at_pos(bits, 5));
    EXPECT_FALSE(bs64::is_set_at_pos(bits, 4));
}

TEST(Bitset64Test, Reset) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 5);
    bits = bs64::clear_bit_at_pos(bits, 5);
    EXPECT_FALSE(bs64::is_set_at_pos(bits, 5));
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

TEST(Bitset64Test, FindFirstZero) {
    bitset64 bits = 0ULL;
    EXPECT_EQ(bs64::find_pos_first_set_bit(bits), 64); // No bits set
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

TEST(Bitset64Test, SmallestRepresentation) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 1);
    bits = bs64::set_bit_at_pos(bits, 2);
    // bits = 0110 (bits 0-3)
    // Rotations: 0110, 0011, 1001, 1100
    // Smallest: 0011
    
    bitset64 smallest = bs64::smallest_representation(bits, 4);
    EXPECT_EQ(bs64::count_set_bits(smallest), 2);
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

TEST(Bitset64Test, NextBitsetWithSamePopcount) {
    bitset64 bits = 0ULL;
    bits = bs64::set_bit_at_pos(bits, 0);
    bits = bs64::set_bit_at_pos(bits, 1);
    bits = bs64::set_bit_at_pos(bits, 2);
    // bits = 00000111
    
    bitset64 next = bs64::next_bitset_with_same_popcount(bits);
    // Next should be 00001011 (bit 2 moves to bit 3)
    EXPECT_EQ(bs64::count_set_bits(next), 3);
    EXPECT_GT(next, bits);
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

// Test Hash Function
TEST(Bitset64Test, Hash) {
    bitset64 bits1 = 0ULL;
    bits1 = bs64::set_bit_at_pos(bits1, 1);
    bits1 = bs64::set_bit_at_pos(bits1, 3);
    
    bitset64 bits2 = 0ULL;
    bits2 = bs64::set_bit_at_pos(bits2, 1);
    bits2 = bs64::set_bit_at_pos(bits2, 3);
    
    std::size_t h1 = bs64::hash(bits1);
    std::size_t h2 = bs64::hash(bits2);
    EXPECT_EQ(h1, h2);
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
    EXPECT_EQ(bs64::find_pos_first_set_bit(bits), 64);
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

