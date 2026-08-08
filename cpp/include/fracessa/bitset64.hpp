// bitset64.hpp
#pragma once

#include <cstdint>
#include <string>
#include <cstddef>

/*
 * A support is the set of pure strategies played with positive probability.
 * Bit i is one exactly when strategy i belongs to that support. This makes set
 * difference, subset tests, and cyclic shifts ordinary integer operations.
 *
 * The complete search may visit 2^n supports, so these small helpers are called
 * millions or billions of times. They are deliberately inline and mostly
 * unchecked; callers must respect the stated dimension and nonzero preconditions.
 */

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

// Caller precondition: x != 0. Deliberately unchecked in this hot primitive.
inline size_t ctz64(uint64_t x) noexcept {
  #ifdef _MSC_VER
    unsigned long index;
    _BitScanForward64(&index, x);
    return static_cast<size_t>(index);
  #else
    return static_cast<size_t>(__builtin_ctzll(x));
  #endif
}

// One word stores every support for dimensions 1 through 64. Complete numeric
// enumeration with a uint64_t one-past-end sentinel still requires n < 64.
typedef uint64_t bitset64;

// Namespace for bitset64 operations
namespace bs64 {

constexpr size_t kMaxBitsetDimension = 64;

// Caller precondition: n < 64 for complete support enumeration.
inline bitset64 two_to_the_power_of(size_t n) noexcept {
  return 1ULL << n;
}

inline bitset64 set_bit_at_pos(bitset64 bits, size_t pos) noexcept {
  return bits | (1ULL << pos);
}

inline bitset64 set_all_n_bits(size_t n) noexcept {
  return n == 64 ? ~bitset64{0} : (bitset64{1} << n) - 1;
}

// Shift every strategy index down by one modulo n. Circular-symmetric games
// use this to obtain the other supports in the same rotational orbit.
inline bitset64 rot_one_right(bitset64 bits, size_t n) noexcept {
  bitset64 mask = set_all_n_bits(n);
  bitset64 low = bits & mask;
  bitset64 lo = low << (n - 1);
  bitset64 hi = low >> 1;
  return (hi | lo) & mask;
}

// Rotate the lowest n bits left by shift positions. Caller guarantees
// 1 <= n <= 64 and shift < n.
inline bitset64 rot_left(bitset64 bits, size_t shift, size_t n) noexcept {
  const bitset64 mask = set_all_n_bits(n);
  bits &= mask;
  if (shift == 0) return bits;
  return ((bits << shift) | (bits >> (n - shift))) & mask;
}

inline bitset64 reverse_bits(bitset64 bits) noexcept {
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
  // GCC does not recognize the portable mask sequence as ARM64's single rbit instruction.
  __asm__("rbit %0, %1" : "=r"(bits) : "r"(bits));
  return bits;
#else
  bits = ((bits >> 1) & 0x5555555555555555ULL) | ((bits & 0x5555555555555555ULL) << 1);
  bits = ((bits >> 2) & 0x3333333333333333ULL) | ((bits & 0x3333333333333333ULL) << 2);
  bits = ((bits >> 4) & 0x0f0f0f0f0f0f0f0fULL) | ((bits & 0x0f0f0f0f0f0f0f0fULL) << 4);
  bits = ((bits >> 8) & 0x00ff00ff00ff00ffULL) | ((bits & 0x00ff00ff00ff00ffULL) << 8);
  bits = ((bits >> 16) & 0x0000ffff0000ffffULL) | ((bits & 0x0000ffff0000ffffULL) << 16);
  return (bits >> 32) | (bits << 32);
#endif
}

// Mirror strategy i to n-1-i. Rotations of this result cover every reflection
// axis of a circular support.
inline bitset64 reflect(bitset64 bits, size_t n) noexcept {
  if (n == 0) return 0;
  return reverse_bits(bits) >> (64 - n);
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

// Return the next set position strictly after pos, or 64 if none remains.
inline size_t find_pos_next_set_bit(bitset64 bits, size_t pos) noexcept {
  size_t p = pos + 1;
  if (p >= 64) return static_cast<size_t>(64);
  uint64_t w = bits & (~0ULL << p);  // zero below p
  if (w) return ctz64(w);
  return static_cast<size_t>(64);  // no more bits found
}

// Set inclusion: bits is a subset of o exactly when it has no bit outside o.
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

inline bitset64 lowest_set_bit_as_bit(bitset64 bits) noexcept {
  // Unsigned subtraction wraps: x & -x isolates the lowest 1 bit and maps 0 to 0.
  return bits & (0ULL - bits);
}

// Gosper's algorithm: return the next larger bit pattern with the same number
// of set bits. Caller guarantees bits != 0 and that another pattern remains in the active width.
inline bitset64 next_same_cardinality(bitset64 bits) noexcept {
  const bitset64 lowest = lowest_set_bit_as_bit(bits);
  const bitset64 ripple = bits + lowest;
  return ripple | (((ripple ^ bits) >> 2) / lowest);
}

// Keep members whose original strategy index is below pos. Counting them maps
// an original strategy index to its row in a compact support-indexed matrix.
inline bitset64 bits_before_pos(bitset64 bits, size_t pos) noexcept {
  return bits & ((1ULL << pos) - 1);
}

// Canonical representative under cyclic rotations:
// for circular symmetry, all n rotations describe the same support orbit.
// We keep only the lexicographically smallest bit pattern in that orbit.
inline bool is_smallest_representation(bitset64 bits, size_t n) noexcept {
  uint64_t mask = set_all_n_bits(n);
  uint64_t original = bits & mask;
  uint64_t current = original;
  size_t shift_left = n - 1;
  
  for (size_t i = 1; i < n; i++) {
    current = ((current >> 1) | (current << shift_left)) & mask;
    if (current < original) {
      return false;
    }
  }
  // All rotations are >= original, so it's canonical
  return true;
}

// Convert to bitstring representation (MSB first, like std::bitset::to_string())
// Only outputs bits 0 to dimension-1 (rightmost dimension bits)
// Example: dimension=5, bits 0,3,4 set -> "11001"
inline std::string to_bitstring(bitset64 bits, size_t dimension) {
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
inline std::string to_string(bitset64 bits) {
  return std::to_string(bits);
}

// Write the set positions in increasing order to caller-owned stack storage.
// ctz finds the current lowest bit; bits &= bits - 1 removes that bit without
// scanning the other 64 positions. The returned count is the support size.
inline size_t extract_set_indices(bitset64 bits, size_t dimension, uint8_t (&indices)[kMaxBitsetDimension]) noexcept
{
  if (dimension < kMaxBitsetDimension) {
    bits &= (1ULL << dimension) - 1ULL;
  }

  size_t count = 0;
  while (bits) {
    const size_t pos = ctz64(bits);
    indices[count++] = static_cast<uint8_t>(pos);
    bits &= bits - 1;
  }
  return count;
}

} // namespace bs64
