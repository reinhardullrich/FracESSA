// bitset64.hpp
#pragma once

#include <cstdint>
#include <string>
#include <cstddef>

// Platform-specific intrinsics
#ifdef _MSC_VER
  #include <intrin.h>
#endif

// Portable popcount wrapper
inline size_t popcount64(uint64_t x) noexcept {
  #ifdef _MSC_VER
    return static_cast<size_t>(_mm_popcnt_u64(x));
  #else
    return static_cast<size_t>(__builtin_popcountll(x));
  #endif
}

// Portable count trailing zeros wrapper
inline size_t ctz64(uint64_t x) noexcept {
  #ifdef _MSC_VER
    unsigned long index;
    if (_BitScanForward64(&index, x)) {
      return static_cast<size_t>(index);
    }
    return 64; // undefined behavior case, but we check for 0 before calling
  #else
    return static_cast<size_t>(__builtin_ctzll(x));
  #endif
}

/// Ultra-optimized bitset for n <= 64.
/// - stores only uint64_t bits (8 bytes, same as uint64_t)
/// - nbits must be provided by caller when needed for masking
/// - all operations are inlined and branch-light
/// - NO bounds checks for maximum performance

// Type alias: bitset64 is just uint64_t
typedef uint64_t bitset64;

// Namespace for bitset64 operations
namespace bs64 {

inline bitset64 two_to_the_power_of(size_t n) noexcept {
  return 1ULL << n;
}

inline bitset64 set_bit_at_pos(bitset64 bits, size_t pos) noexcept {
  return bits | (1ULL << pos);
}

inline bitset64 clear_bit_at_pos(bitset64 bits, size_t pos) noexcept {
  return bits & ~(1ULL << pos);
}

inline bitset64 set_all_n_bits(size_t n) noexcept {
  return (1ULL << n) - 1ULL; //careful with n == 0
}

// circular rotate right by exactly 1 bit for a bitmask of size nbits
inline bitset64 rot_one_right(bitset64 bits, size_t n) noexcept {
  bitset64 mask = set_all_n_bits(n);
  bitset64 low = bits & mask;
  bitset64 lo = low << (n - 1);
  bitset64 hi = low >> 1;
  return (hi | lo) & mask;
}

inline bool is_set_at_pos(bitset64 bits, size_t pos) noexcept {
  return (bits >> pos) & 1ULL;
}

inline size_t count_set_bits(bitset64 bits) noexcept {
  return popcount64(bits);
}

inline size_t find_pos_first_set_bit(bitset64 bits) noexcept {
  return ctz64(bits);
}

// find next bit after pos
inline size_t find_pos_next_set_bit(bitset64 bits, size_t pos) noexcept {
  size_t p = pos + 1;
  uint64_t w = bits & (~0ULL << p);  // zero below p
  if (w) return ctz64(w);
  return static_cast<size_t>(64);  // no more bits found
}

// this ⊆ o <=> (this & ~o) == 0
inline bool is_subset_of(bitset64 bits, bitset64 o) noexcept {
  return (bits & ~o) == 0ULL;
}

// Return this \ o (set difference)
inline bitset64 subtract(bitset64 bits, bitset64 o) noexcept {
  return bits & ~o;
}

// Get a single-bit bitset64 with only the bit at position pos set
inline bitset64 single_bit_at_pos(size_t pos) noexcept {
  return 1ULL << pos;
}

// Get the lowest set bit as a bitset64 (only that bit is set, all other bits are 0)
inline bitset64 lowest_set_bit_as_bit(bitset64 bits) noexcept {
  size_t pos = find_pos_first_set_bit(bits);
  return single_bit_at_pos(pos);
}

// Get all set bits in bits that are before position pos
inline bitset64 bits_before_pos(bitset64 bits, size_t pos) noexcept {
  return bits & ((1ULL << pos) - 1);
}

//careful! no check for the biggest element here! done outside!
inline bitset64 next_bitset_with_same_popcount(bitset64 bits) noexcept {
  uint64_t t = bits | (bits - 1);
  return (t + 1) | (((~t & -~t) - 1) >> (ctz64(bits) + 1));
}

// Check if this bitset is in its smallest representation (canonical form)
// Uses early exit: returns false immediately if any rotation is smaller than original
// Optimized version: reduces masking operations for better performance
inline bool is_smallest_representation(bitset64 bits, size_t n) noexcept {
  uint64_t mask = set_all_n_bits(n);
  uint64_t original = bits & mask;
  uint64_t current = original;
  size_t shift_left = n - 1;
  
  for (size_t i = 1; i < n; i++) {
    current = ((current >> 1) | (current << shift_left)) & mask;
    if (current < original) {
      return false; //there exists a smaller representation, ie. this one cannot be canonical!
    }
  }
  // All rotations are >= original, so it's canonical
  return true;
}

// Convert to bitstring representation (MSB first, like std::bitset::to_string())
// Only outputs bits 0 to dimension-1 (rightmost dimension bits)
// Example: dimension=5, bits 0,3,4 set -> "10011"
inline std::string to_bitstring(bitset64 bits, size_t dimension) noexcept {
  if (dimension == 0) return "";
  
  std::string result;
  result.reserve(dimension);
  
  // Output from highest bit (dimension-1) to lowest bit (0)
  for (int i = static_cast<int>(dimension) - 1; i >= 0; i--) {
    result += (is_set_at_pos(bits, static_cast<size_t>(i)) ? '1' : '0');
  }
  
  return result;
}

// Convert to string representation (decimal representation of uint64)
inline std::string to_string(bitset64 bits) noexcept {
  return std::to_string(bits);
}

// portable hash (fast)
inline std::size_t hash(bitset64 bits) noexcept {
  // 64-bit mix (FNV-like)
  uint64_t x = bits;
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return static_cast<std::size_t>(x);
}

} // namespace bs64