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
    EXPECT_EQ(bits, bitset{1} << 5);
}

TEST(BitsetTest, SetAll) {
    bitset bits = 0ULL;
    bits = set_all_n_bits(5);
    EXPECT_EQ(bits, 0b11111u);
    EXPECT_EQ(count_set_bits(bits), 5);
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
    EXPECT_EQ(bits, 0b10010u);
}

TEST(BitsetTest, RotOneRightMasksBitsOutsideDimension) {
    // Bit 6 is outside n=5 and must be masked out.
    bitset bits = 0ULL;
    bits = set_bit_at_pos(bits, 0);
    bits = set_bit_at_pos(bits, 6);

    bits = rot_one_right(bits, 5);

    EXPECT_EQ(bits, bitset{1} << 4);
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
    EXPECT_EQ(bits, 0xffu);
    EXPECT_EQ(count_set_bits(bits), 8);
}

TEST(BitsetTest, SetAllZero) {
    bitset bits = 0ULL;
    bits = set_all_n_bits(0);
    EXPECT_EQ(count_set_bits(bits), 0);
}
