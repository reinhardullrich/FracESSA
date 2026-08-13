#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <fracessa/bitset.hpp>

namespace fracessa::support {

/*
 * Independent test oracle for circular support generation. The FKM recursion creates one mask for all rotations; the reflection
 * check then removes the duplicate mirror image. Whenever one representative is ready, the generator immediately calls:
 *
 *     callback(support, support_size)
 *
 * The callback order is cardinality first; within one cardinality, the emitted representative masks increase numerically.
 *
 * add_forbidden_orbit() stores every rotation and reflection of an exact equilibrium support. This lets the recursion prune
 * equivalent forbidden subsets in every orientation while it generates later support sizes. Its return value is the
 * number of distinct masks in that dihedral orbit: the candidate multiplier.
 */
class ReferenceCircularSupportGenerator {
private:
    size_t dimension_;
    size_t target_cardinality_ = 0;
    // The FKM recursion uses this one-indexed array as its working binary word.
    std::array<uint8_t, kMaxBitsetDimension + 1> word_{};
    // A rule is checked when its lowest set bit is added, because that is the first moment the partial support can
    // contain the complete rule.
    std::array<std::vector<bitset>, kMaxBitsetDimension> forbidden_by_lowest_;
    // Rules found by the callback stay pending until the next support size.
    std::vector<bitset> pending_forbidden_;
    bool emitted_ = false;

    inline void activate_pending() {
        for (const bitset support : pending_forbidden_)
            forbidden_by_lowest_[ctz64(support)].push_back(support);
        pending_forbidden_.clear();
    }

    inline bool completes_forbidden(bitset partial, size_t new_lowest_bit) const noexcept {
        for (const bitset forbidden : forbidden_by_lowest_[new_lowest_bit]) {
            if (is_subset_of(forbidden, partial))
                return true;
        }
        return false;
    }

    inline bitset smallest_rotation(bitset support) const noexcept {
        bitset smallest = support;
        for (size_t i = 1; i < dimension_; ++i) {
            support = rot_one_right(support, dimension_);
            smallest = std::min(smallest, support);
        }
        return smallest;
    }

    template<class Callback>
    inline void emit_necklace(bitset support, Callback& callback) {
        if (support <= smallest_rotation(reflect(support, dimension_))) {
            emitted_ = true;
            callback(support, target_cardinality_);
        }
    }

    template<class Callback>
    inline void generate_necklaces(size_t position, size_t period, size_t ones, bitset partial, Callback& callback) {
        if (position > dimension_) {
            if (dimension_ % period == 0 && ones == target_cardinality_)
                emit_necklace(partial, callback);
            return;
        }
        if (ones > target_cardinality_ || ones + dimension_ - position + 1 < target_cardinality_)
            return;

        const size_t bit = dimension_ - position;
        word_[position] = word_[position - period];
        if (word_[position] == 0) {
            generate_necklaces(position + 1, period, ones, partial, callback);
        } else {
            const bitset with_bit = set_bit_at_pos(partial, bit);
            if (!completes_forbidden(with_bit, bit)) {
                generate_necklaces(position + 1, period, ones + 1, with_bit, callback);
            }
        }

        if (word_[position - period] == 0) {
            word_[position] = 1;
            const bitset with_bit = set_bit_at_pos(partial, bit);
            if (!completes_forbidden(with_bit, bit)) {
                generate_necklaces(position + 1, position, ones + 1, with_bit, callback);
            }
        }
    }

public:
    explicit ReferenceCircularSupportGenerator(size_t dimension) noexcept : dimension_(dimension) {}

    // Callback signature: void(bitset support, size_t support_size).
    // It is called synchronously once for each generated representative.
    template<class Callback>
    inline void generate(Callback&& callback) {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            word_.fill(0);
            emitted_ = false;
            generate_necklaces(1, 1, 0, 0, callback);
            // If no support of this size survives, no larger support can survive either: every larger support contains
            // one of this size.
            if (!emitted_)
                return;
        }
    }

    inline size_t add_forbidden_orbit(bitset support) {
        const bitset reflected = reflect(support, dimension_);
        bool reflection_is_rotation = false;
        size_t multiplier = 0;
        bitset current = support;
        do {
            reflection_is_rotation = reflection_is_rotation || current == reflected;
            pending_forbidden_.push_back(current);
            ++multiplier;
            current = rot_one_right(current, dimension_);
        } while (current != support);

        if (!reflection_is_rotation) {
            current = reflected;
            do {
                pending_forbidden_.push_back(current);
                ++multiplier;
                current = rot_one_right(current, dimension_);
            } while (current != reflected);
        }
        return multiplier;
    }
};

} // namespace fracessa::support
