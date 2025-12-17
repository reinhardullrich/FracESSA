// supports.hpp
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fracessa/bitset64.hpp>

// Force-inline hint (same as bitset64.hpp)
#if defined(_MSC_VER)
#  define FORCE_INLINE __forceinline
#else
#  define FORCE_INLINE __attribute__((always_inline)) inline
#endif

/// Compute binomial coefficient C(n,k) = n!/(k!(n-k)!)
/// Returns uint64_t - safe since n <= 64, no overflow possible
/// No safety checks needed per user requirement
inline uint64_t binomial_coefficient(uint64_t n, uint64_t k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k; // Use symmetry
    
    uint64_t result = 1;
    for (uint64_t i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

/// High-performance Supports class for managing support sets
/// All hot-path methods are FORCE_INLINE for maximum performance
class Supports {
private:
    std::vector<std::vector<bitset64>> supports_;
    size_t dimension_;
    bool is_cs_;
    
public:
    /// Constructor - stores parameters, does not initialize supports
    inline Supports(size_t dimension, bool is_cs) noexcept
        : supports_(dimension)
        , dimension_(dimension)
        , is_cs_(is_cs)
    {}
    
    /// Initialize all supports - must be called after construction
    inline void initialize() {
        // Reserve space for each support size using binomial coefficients and set coprime flags
        std::vector<bool> is_coprime(dimension_);
        for (size_t i = 0; i < dimension_; ++i) {
            supports_[i].reserve(binomial_coefficient(dimension_, i + 1));            
            if (is_cs_) 
                    is_coprime[i] = (std::gcd(i+1, dimension_) == 1);
        }
        
        // Populate supports based on is_cs_ flag
        if (is_cs_) {
            for (uint64_t support = 1ULL; support < bs64::two_to_the_power_of(dimension_); ++support) {
                size_t current_index = bs64::count_set_bits(support) - 1;
                if (is_coprime[current_index]) {
                    // Only add if it's the smallest representation (canonical form)
                    if (bs64::is_smallest_representation(support, dimension_)) {
                        supports_[current_index].push_back(support);
                    }
                } else {
                    supports_[current_index].push_back(support);
                }
            }
        } else {
            for (uint64_t bits = 1ULL; bits < bs64::two_to_the_power_of(dimension_); ++bits) {
                bitset64 support = bits;
                supports_[bs64::count_set_bits(support) - 1].push_back(support);
            }
        }
    }
    
    /// Get const reference to supports for a given support size (1-indexed)
    /// CRITICAL hot path - must be FORCE_INLINE
    inline const std::vector<bitset64>& get_supports(size_t support_size) const noexcept {
        return supports_[support_size-1];
    }
    
    /// Remove all supersets of the given subset, starting from from_size
    /// Hot path - FORCE_INLINE for maximum performance
    inline void remove_supersets(const bitset64& subset, uint64_t support_size = 0) noexcept {
        if (support_size == 0) {
            support_size = bs64::count_set_bits(subset);
        }
        for (size_t i = support_size; i < dimension_; ++i) { //index support_size means erase from support_size+1 on!!!
            auto& vec = supports_[i];
            // Binary search: find first element where x > subset (only for dimension >= 10)
            // this is an educated guess. for small dimensions the overhead of the binary search is bigger than the time saved
            auto start_it = (dimension_ >= 10) 
                ? std::upper_bound(vec.begin(), vec.end(), subset)
                : vec.begin();
            // Only check elements from start_it onwards
            if (start_it != vec.end()) {
                vec.erase(
                std::remove_if(
                        start_it,
                        vec.end(),
                    [=](const bitset64& x) { return bs64::is_subset_of(subset, x); }
                ),
                    vec.end()
            );
            }
        }
    }
    
};

