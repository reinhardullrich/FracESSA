// Experimental circular support generator with dihedral-orbit pruning.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include <fracessa/bitset64.hpp>

struct SupportOrbit {
    // Number of different supports obtained by rotation alone.
    size_t rotation_count = 1;
    // True when mirroring produces another, disjoint set of rotations.
    bool has_reflected_class = false;

    size_t member_count() const noexcept {
        return rotation_count * (has_reflected_class ? 2 : 1);
    }
};

// Mirror positions 0,...,n-1 by mapping i to n-1-i. Rotations of this result
// cover every possible reflection axis.
inline bitset64 reflect_support(bitset64 support, size_t dimension) noexcept {
    bitset64 reflected = 0;
    for (size_t i = 0; i < dimension; ++i)
        reflected |= ((support >> i) & 1ULL) << (dimension - 1 - i);
    return reflected;
}

inline uint64_t binomial_coefficient(uint64_t n, uint64_t k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;

    uint64_t result = 1;
    for (uint64_t i = 0; i < k; ++i)
        result = result * (n - i) / (i + 1);
    return result;
}

/*
 * In a circular game, rotating or mirroring a support gives an equivalent
 * mathematical case. All supports related this way form a "bracelet". We
 * generate only one support from each bracelet instead of generating all 2^n
 * masks and discarding equivalent ones afterward.
 *
 * Once a support has produced an exact equilibrium, no strict superset can be
 * an ESS. For circular games we remember that support and all its rotations
 * and reflections, then skip future generated supports containing one of them.
 * Non-circular games retain the production implementation, which stores all
 * supports and physically erases such supersets.
 */
class Supports {
private:
    std::vector<std::vector<bitset64>> supports_;
    std::vector<bitset64> generated_;
    // Store each earlier candidate support under its lowest set bit. A future
    // support can contain it only if that bit is also present, so unrelated
    // buckets never need to be checked.
    std::vector<std::vector<bitset64>> forbidden_by_anchor_;
    std::vector<uint8_t> word_;
    size_t dimension_;
    size_t target_size_ = 0;
    bool is_cs_;
    bool exhausted_ = false;

    bitset64 smallest_rotation(bitset64 support) const noexcept {
        bitset64 smallest = support;
        for (size_t i = 1; i < dimension_; ++i) {
            support = bs64::rot_one_right(support, dimension_);
            smallest = std::min(smallest, support);
        }
        return smallest;
    }

    bool is_pruned(bitset64 support) const noexcept {
        bitset64 anchors = support;
        while (anchors != 0) {
            const size_t anchor = ctz64(anchors);
            for (const bitset64 subset : forbidden_by_anchor_[anchor]) {
                if (bs64::is_subset_of(subset, support))
                    return true;
            }
            anchors &= anchors - 1;
        }
        return false;
    }

    void emit_necklace() {
        bitset64 support = 0;
        for (size_t i = 1; i <= dimension_; ++i) {
            if (word_[i] != 0)
                support |= 1ULL << (dimension_ - i);
        }

        // FKM already gives the smallest rotation. Mirroring can give a second
        // rotation family; keep whichever family has the smaller representative.
        const bitset64 reflected = reflect_support(support, dimension_);
        if (support <= smallest_rotation(reflected) && !is_pruned(support))
            generated_.push_back(support);
    }

    // Generate binary necklaces containing exactly target_size_ set bits. A
    // necklace represents all rotations of one support, so this recursion does
    // not visit every one of the C(n,k) individual masks.
    void generate_necklaces(size_t position, size_t period, size_t ones) {
        if (position > dimension_) {
            if (dimension_ % period == 0 && ones == target_size_)
                emit_necklace();
            return;
        }
        if (ones > target_size_
            || ones + dimension_ - position + 1 < target_size_)
            return;

        word_[position] = word_[position - period];
        generate_necklaces(position + 1, period, ones + word_[position]);

        if (word_[position - period] == 0) {
            word_[position] = 1;
            generate_necklaces(position + 1, position, ones + 1);
        }
    }

    void add_forbidden(bitset64 support) {
        forbidden_by_anchor_[ctz64(support)].push_back(support);
    }

public:
    inline Supports(size_t dimension, bool is_cs) noexcept
        : supports_(dimension)
        , forbidden_by_anchor_(dimension)
        , word_(dimension + 1, 0)
        , dimension_(dimension)
        , is_cs_(is_cs)
    {}

    inline void initialize() {
        generated_.clear();
        for (auto& bucket : supports_)
            bucket.clear();
        for (auto& bucket : forbidden_by_anchor_)
            bucket.clear();
        exhausted_ = false;

        if (is_cs_)
            return;

        // Non-circular control path: preserve production's eager mask order.
        for (size_t i = 0; i < dimension_; ++i)
            supports_[i].reserve(binomial_coefficient(dimension_, i + 1));
        for (bitset64 support = 1;
             support < bs64::two_to_the_power_of(dimension_); ++support) {
            supports_[bs64::count_set_bits(support) - 1].push_back(support);
        }
    }

    inline const std::vector<bitset64>& get_supports(size_t support_size) {
        if (!is_cs_)
            return supports_[support_size - 1];

        generated_.clear();
        if (exhausted_)
            return generated_;

        target_size_ = support_size;
        std::fill(word_.begin(), word_.end(), 0);
        generate_necklaces(1, 1, 0);

        // If every size-k support was excluded, every larger support contains
        // an excluded size-k support. The search can therefore stop completely.
        exhausted_ = generated_.empty();
        return generated_;
    }

    inline SupportOrbit remove_supersets(
        const bitset64& subset, uint64_t support_size = 0) {
        if (support_size == 0)
            support_size = bs64::count_set_bits(subset);

        if (!is_cs_) {
            // Here future supports already exist in memory, so remove every
            // strict superset immediately, as the production code does.
            for (size_t i = support_size; i < dimension_; ++i) {
                auto& bucket = supports_[i];
                const auto begin = dimension_ >= 10
                    ? std::upper_bound(bucket.begin(), bucket.end(), subset)
                    : bucket.begin();
                bucket.erase(
                    std::remove_if(
                        begin, bucket.end(),
                        [=](bitset64 support) {
                            return bs64::is_subset_of(subset, support);
                        }),
                    bucket.end());
            }
            return {};
        }

        /*
         * Circular supports are generated later, so there are no stored
         * supersets to erase. Instead, record every rotated and mirrored copy
         * of this exact equilibrium support. Future generation uses these as
         * subset tests and skips the corresponding supersets.
         *
         * This same traversal tells the caller how many equivalent candidates
         * exist. Returning that information avoids walking the orbit a second
         * time when candidate rows and the ESS count are reconstructed.
         */
        SupportOrbit orbit;
        orbit.rotation_count = 0;
        const bitset64 reflected = reflect_support(subset, dimension_);
        bool reflection_is_rotation = false;

        bitset64 current = subset;
        do {
            reflection_is_rotation = reflection_is_rotation || current == reflected;
            add_forbidden(current);
            ++orbit.rotation_count;
            current = bs64::rot_one_right(current, dimension_);
        } while (current != subset);

        orbit.has_reflected_class = !reflection_is_rotation;
        if (orbit.has_reflected_class) {
            current = reflected;
            do {
                add_forbidden(current);
                current = bs64::rot_one_right(current, dimension_);
            } while (current != reflected);
        }
        return orbit;
    }
};
