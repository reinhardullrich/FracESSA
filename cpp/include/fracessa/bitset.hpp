// bitset.hpp
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

namespace fracessa::support {

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

// One word stores every support for dimensions 1 through 64; dimension 64 uses all bits. Larger dimensions use bitset_multiword,
// while this type remains the allocation-free hot path for small games.
using bitset = uint64_t;

constexpr size_t kMaxBitsetDimension = 64;

inline bitset set_bit_at_pos(bitset bits, size_t pos) noexcept {
  return bits | (1ULL << pos);
}

inline bitset set_all_n_bits(size_t n) noexcept {
  return n == 64 ? ~bitset{0} : (bitset{1} << n) - 1;
}

// Shift every strategy index down by one modulo n. Circular-symmetric games
// use this to obtain the other supports in the same rotational orbit.
inline bitset rot_one_right(bitset bits, size_t n) noexcept {
  bitset mask = set_all_n_bits(n);
  bitset low = bits & mask;
  bitset lo = low << (n - 1);
  bitset hi = low >> 1;
  return (hi | lo) & mask;
}

inline bitset reverse_bits(bitset bits) noexcept {
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
inline bitset reflect(bitset bits, size_t n) noexcept {
  if (n == 0) return 0;
  return reverse_bits(bits) >> (64 - n);
}

inline size_t count_set_bits(bitset bits) noexcept {
  return popcount64(bits);
}

inline size_t find_pos_first_set_bit(bitset bits) noexcept {
  return ctz64(bits);
}

// Set inclusion: bits is a subset of o exactly when it has no bit outside o.
inline bool is_subset_of(bitset bits, bitset o) noexcept {
  return (bits & ~o) == 0ULL;
}

// Return this \ o (set difference)
inline bitset subtract(bitset bits, bitset o) noexcept {
  return bits & ~o;
}

// Convert to string representation (decimal representation of uint64)
inline std::string to_string(bitset bits) {
  return std::to_string(bits);
}

// Format the set positions in increasing order for diagnostic output.
inline std::string to_index_set(bitset bits) {
  std::string result = "{";
  bool first = true;
  while (bits != 0) {
    if (!first) result += ", ";
    result += std::to_string(ctz64(bits));
    first = false;
    bits &= bits - 1;
  }
  result += '}';
  return result;
}

// Write the set positions in increasing order to caller-owned stack storage.
// ctz finds the current lowest bit; bits &= bits - 1 removes that bit without
// scanning the other 64 positions. The returned count is the support size.
inline size_t extract_set_indices(bitset bits, size_t dimension, uint8_t (&indices)[kMaxBitsetDimension]) noexcept
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

} // namespace fracessa::support
