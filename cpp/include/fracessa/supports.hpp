// supports.hpp
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fracessa/bitset64.hpp>

/*
 * Stores the remaining support sets, grouped by cardinality.
 * `supports_[k-1]` contains the masks of size k. The search visits these groups
 * from small to large because an exact equilibrium on support I rules out every
 * larger ESS support containing I.
 *
 * The pruning follows Bomze (1992), Proposition 2.1. If x is an equilibrium
 * and p is an ESS with I(x) contained in J(p), then x = p. In particular, an ESS
 * cannot have a support that strictly contains I(x), because I(p) is contained
 * in J(p). This is why pruning happens after every exact candidate, not only
 * after a candidate has passed the stability test.
 */

// Compute C(n,k) as a vector-capacity hint. Search correctness does not depend
// on this reservation estimate; normal callers satisfy 0 <= k <= n.
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

// The constructor creates empty cardinality buckets. Enumeration is delayed so
// --fullsupport can test one mask without first materializing all 2^n supports.
class Supports {
private:
    std::vector<std::vector<bitset64>> supports_;
    size_t dimension_;
    bool is_cs_;
    
public:
    inline Supports(size_t dimension, bool is_cs) noexcept
        : supports_(dimension)
        , dimension_(dimension)
        , is_cs_(is_cs)
    {}
    
    /*
     * Enumerate every nonempty support in increasing integer-mask order.
     *
     * In a circular-symmetric game, rotated masks describe equivalent cases. If
     * gcd(k,n)=1, a size-k mask cannot have a shorter rotational period, so its
     * orbit contains exactly n distinct masks. We store only the smallest mask
     * and reconstruct its other n-1 rotations after solving it. For noncoprime
     * k, shorter orbits are possible, so this simple reduction is not used.
     */
    inline void initialize() {
        // Reserve once before enumeration; these vectors can otherwise reallocate
        // repeatedly while collecting a large support family.
        std::vector<bool> is_coprime(dimension_);
        for (size_t i = 0; i < dimension_; ++i) {
            supports_[i].reserve(binomial_coefficient(dimension_, i + 1));            
            if (is_cs_) 
                    is_coprime[i] = (std::gcd(i+1, dimension_) == 1);
        }
        
        if (is_cs_) {
            for (uint64_t support = 1ULL; support < bs64::two_to_the_power_of(dimension_); ++support) {
                size_t current_index = bs64::count_set_bits(support) - 1;
                if (is_coprime[current_index]) {
                    // Keep one representative of this full-size rotation orbit.
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
    
    // support_size is mathematical cardinality, hence the minus-one bucket index.
    inline const std::vector<bitset64>& get_supports(size_t support_size) const noexcept {
        return supports_[support_size-1];
    }
    
    // Remove strict supersets from every larger cardinality bucket. `subset` is
    // known to support an exact equilibrium, which is the pruning precondition
    // explained above.
    inline void remove_supersets(const bitset64& subset, uint64_t support_size = 0) noexcept {
        if (support_size == 0) {
            support_size = bs64::count_set_bits(subset);
        }
        // Bucket i contains sets of size i+1, so starting at support_size skips
        // the candidate's own cardinality and begins with strict supersets.
        for (size_t i = support_size; i < dimension_; ++i) {
            auto& vec = supports_[i];
            // Masks were inserted in numeric order, and every superset is larger
            // than `subset`. The n=10 cutoff is a speed heuristic: above it the
            // code skips the impossible prefix, below it the scan starts directly.
            auto start_it = (dimension_ >= 10) 
                ? std::upper_bound(vec.begin(), vec.end(), subset)
                : vec.begin();
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
