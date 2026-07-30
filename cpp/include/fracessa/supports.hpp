// supports.hpp
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <fracessa/bitset64.hpp>

// Kept here because the exact copositivity code also uses it.
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
 * Generate every nonempty support, one mask at a time, without storing a list. generate() visits support sizes 1 through
 * dimension. For each size, a binary depth-first search decides which bits belong to the support. Whenever it completes
 * one support, it immediately calls the supplied callback with:
 *
 *     callback(support, support_size)
 *
 * The callback order is part of the contract: support_size increases first; within one support size, the uint64_t
 * support masks increase numerically.
 *
 * add_forbidden() records an exact equilibrium support. Before generating the next size, that support becomes a pruning
 * rule: every recursive branch that already contains it is skipped because all masks below that branch are also
 * supersets of it.
 */
class NonCircularSupportGenerator {
private:
    size_t dimension_;
    // A rule is checked when its lowest set bit is added, because that is the first moment the partial support can
    // contain the complete rule.
    std::array<std::vector<bitset64>, bs64::kMaxBitsetDimension> forbidden_by_lowest_;
    // Rules found by the callback stay pending until the next support size.
    std::vector<bitset64> pending_forbidden_;
    size_t target_cardinality_ = 0;
    bool emitted_ = false;

    inline void activate_pending() {
        for (const bitset64 support : pending_forbidden_)
            forbidden_by_lowest_[ctz64(support)].push_back(support);
        pending_forbidden_.clear();
    }

    inline bool completes_forbidden(bitset64 partial, size_t new_lowest_bit) const noexcept {
        for (const bitset64 forbidden : forbidden_by_lowest_[new_lowest_bit]) {
            if (bs64::is_subset_of(forbidden, partial))
                return true;
        }
        return false;
    }

    template<class Callback>
    inline void generate_from(size_t bits_remaining, size_t needed, bitset64 partial, Callback& callback) {
        if (needed == 0) {
            emitted_ = true;
            callback(partial, target_cardinality_);
            return;
        }
        if (needed > bits_remaining)
            return;

        const size_t bit = bits_remaining - 1;
        if (needed < bits_remaining)
            generate_from(bit, needed, partial, callback);

        const bitset64 with_bit = bs64::set_bit_at_pos(partial, bit);
        if (!completes_forbidden(with_bit, bit))
            generate_from(bit, needed - 1, with_bit, callback);
    }

public:
    explicit NonCircularSupportGenerator(size_t dimension) noexcept : dimension_(dimension) {}

    // Callback signature: void(bitset64 support, size_t support_size).
    // It is called synchronously once for each generated support.
    template<class Callback>
    inline void generate(Callback&& callback) {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            emitted_ = false;
            generate_from(dimension_, target_cardinality_, 0, callback);
            // If no support of this size survives, no larger support can survive either: every larger support contains
            // one of this size.
            if (!emitted_)
                return;
        }
    }

    inline void add_forbidden(bitset64 support) {
        pending_forbidden_.push_back(support);
    }
};

/*
 * Generate circular supports through the same callback contract as the class above, but emit only one representative
 * for equivalent masks. The FKM recursion creates one mask for all of its rotations; the reflection check then removes
 * the duplicate mirror image. Whenever one representative is ready, the generator immediately calls:
 *
 *     callback(support, support_size)
 *
 * The callback order is the same as above: support_size increases first; within one support size, the emitted
 * representative masks increase numerically.
 *
 * add_forbidden() stores every rotation and reflection of an exact equilibrium support. This lets the recursion prune
 * equivalent forbidden subsets in every orientation while it generates later support sizes. Its return value is the
 * number of distinct masks in that dihedral orbit: the candidate multiplier.
 */
class CircularSupportGenerator {
private:
    size_t dimension_;
    size_t target_cardinality_ = 0;
    // The FKM recursion uses this one-indexed array as its working binary word.
    std::array<uint8_t, bs64::kMaxBitsetDimension + 1> word_{};
    // A rule is checked when its lowest set bit is added, because that is the first moment the partial support can
    // contain the complete rule.
    std::array<std::vector<bitset64>, bs64::kMaxBitsetDimension> forbidden_by_lowest_;
    // Rules found by the callback stay pending until the next support size.
    std::vector<bitset64> pending_forbidden_;
    bool emitted_ = false;

    inline void activate_pending() {
        for (const bitset64 support : pending_forbidden_)
            forbidden_by_lowest_[ctz64(support)].push_back(support);
        pending_forbidden_.clear();
    }

    inline bool completes_forbidden(bitset64 partial, size_t new_lowest_bit) const noexcept {
        for (const bitset64 forbidden : forbidden_by_lowest_[new_lowest_bit]) {
            if (bs64::is_subset_of(forbidden, partial))
                return true;
        }
        return false;
    }

    inline bitset64 smallest_rotation(bitset64 support) const noexcept {
        bitset64 smallest = support;
        for (size_t i = 1; i < dimension_; ++i) {
            support = bs64::rot_one_right(support, dimension_);
            smallest = std::min(smallest, support);
        }
        return smallest;
    }

    template<class Callback>
    inline void emit_necklace(bitset64 support, Callback& callback) {
        if (support <= smallest_rotation(bs64::reflect(support, dimension_))) {
            emitted_ = true;
            callback(support, target_cardinality_);
        }
    }

    template<class Callback>
    inline void generate_necklaces(size_t position, size_t period, size_t ones, bitset64 partial, Callback& callback) {
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
            const bitset64 with_bit = bs64::set_bit_at_pos(partial, bit);
            if (!completes_forbidden(with_bit, bit)) {
                generate_necklaces(position + 1, period, ones + 1, with_bit, callback);
            }
        }

        if (word_[position - period] == 0) {
            word_[position] = 1;
            const bitset64 with_bit = bs64::set_bit_at_pos(partial, bit);
            if (!completes_forbidden(with_bit, bit)) {
                generate_necklaces(position + 1, position, ones + 1, with_bit, callback);
            }
        }
    }

public:
    explicit CircularSupportGenerator(size_t dimension) noexcept : dimension_(dimension) {}

    // Callback signature: void(bitset64 support, size_t support_size).
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

    inline size_t add_forbidden(bitset64 support) {
        const bitset64 reflected = bs64::reflect(support, dimension_);
        bool reflection_is_rotation = false;
        size_t multiplier = 0;
        bitset64 current = support;
        do {
            reflection_is_rotation = reflection_is_rotation || current == reflected;
            pending_forbidden_.push_back(current);
            ++multiplier;
            current = bs64::rot_one_right(current, dimension_);
        } while (current != support);

        if (!reflection_is_rotation) {
            current = reflected;
            do {
                pending_forbidden_.push_back(current);
                ++multiplier;
                current = bs64::rot_one_right(current, dimension_);
            } while (current != reflected);
        }
        return multiplier;
    }
};

/*
 * Experimental test-only alternative to CircularSupportGenerator. Production does not use V2. A later generator
 * comparison test will exercise both versions and decide which implementation to keep once that test framework is ready.
 *
 * V2 stores each forbidden support once. During recursion, two uint64_t masks test all rotational and reflected
 * alignments in parallel instead of storing every mask in the support's dihedral orbit.
 */
class CircularSupportGeneratorV2 {
private:
    size_t dimension_;
    bitset64 dimension_mask_;
    size_t target_cardinality_ = 0;
    std::array<uint8_t, bs64::kMaxBitsetDimension + 1> word_{};
    // Each forbidden support is stored in exactly one cardinality bucket.
    std::array<std::vector<bitset64>, bs64::kMaxBitsetDimension> forbidden_by_cardinality_;
    std::vector<bitset64> pending_forbidden_;
    size_t emitted_multiplier_ = 0;
    bool emitted_ = false;

    inline void activate_pending() {
        for (const bitset64 support : pending_forbidden_)
            forbidden_by_cardinality_[bs64::count_set_bits(support)].push_back(support);
        pending_forbidden_.clear();
    }

    /*
     * Each bit in possible_rotations and possible_reflections represents one circular shift. For every bit in the
     * forbidden support, intersect away shifts whose aligned bit is absent from the partial support. A surviving shift
     * means one rotated or reflected copy is already contained.
     *
     * negated_partial stores position -i modulo dimension for every selected position i. Keeping it beside partial makes
     * each alignment mask one rotate-and-intersect operation rather than an n-step orbit scan.
     */
    inline bool some_orbit_member_is_subset(bitset64 forbidden, bitset64 negated_partial) const noexcept {
        bitset64 possible_rotations = dimension_mask_;
        bitset64 possible_reflections = dimension_mask_;

        while (forbidden) {
            const size_t bit = ctz64(forbidden);
            possible_rotations &= bs64::rot_left(negated_partial, bit, dimension_);
            possible_reflections &= bs64::rot_left(negated_partial, dimension_ - 1 - bit, dimension_);
            if ((possible_rotations | possible_reflections) == 0)
                return false;
            forbidden &= forbidden - 1;
        }
        return true;
    }

    inline bool completes_forbidden(bitset64 negated_partial, size_t partial_cardinality) const noexcept {
        for (size_t cardinality = 1; cardinality <= partial_cardinality; ++cardinality) {
            for (const bitset64 forbidden : forbidden_by_cardinality_[cardinality]) {
                if (some_orbit_member_is_subset(forbidden, negated_partial))
                    return true;
            }
        }
        return false;
    }

    inline bitset64 smallest_rotation(bitset64 support) const noexcept {
        bitset64 smallest = support;
        for (size_t i = 1; i < dimension_; ++i) {
            support = bs64::rot_one_right(support, dimension_);
            smallest = std::min(smallest, support);
        }
        return smallest;
    }

    template<class Callback>
    inline void emit_necklace(bitset64 support, size_t period, Callback& callback) {
        const bitset64 reflected_representative = smallest_rotation(bs64::reflect(support, dimension_));
        if (support <= reflected_representative) {
            emitted_multiplier_ = period * (support == reflected_representative ? 1 : 2);
            emitted_ = true;
            callback(support, target_cardinality_);
        }
    }

    template<class Callback>
    inline void generate_necklaces(size_t position, size_t period, size_t ones, bitset64 partial,
                                   bitset64 negated_partial, Callback& callback) {
        if (position > dimension_) {
            if (dimension_ % period == 0 && ones == target_cardinality_)
                emit_necklace(partial, period, callback);
            return;
        }
        if (ones > target_cardinality_ || ones + dimension_ - position + 1 < target_cardinality_)
            return;

        const size_t bit = dimension_ - position;
        word_[position] = word_[position - period];
        if (word_[position] == 0) {
            generate_necklaces(position + 1, period, ones, partial, negated_partial, callback);
        } else {
            const bitset64 with_bit = bs64::set_bit_at_pos(partial, bit);
            const size_t negated_bit = bit == 0 ? 0 : dimension_ - bit;
            const bitset64 with_negated_bit = bs64::set_bit_at_pos(negated_partial, negated_bit);
            if (!completes_forbidden(with_negated_bit, ones + 1)) {
                generate_necklaces(position + 1, period, ones + 1, with_bit, with_negated_bit, callback);
            }
        }

        if (word_[position - period] == 0) {
            word_[position] = 1;
            const bitset64 with_bit = bs64::set_bit_at_pos(partial, bit);
            const size_t negated_bit = bit == 0 ? 0 : dimension_ - bit;
            const bitset64 with_negated_bit = bs64::set_bit_at_pos(negated_partial, negated_bit);
            if (!completes_forbidden(with_negated_bit, ones + 1)) {
                generate_necklaces(position + 1, position, ones + 1, with_bit, with_negated_bit, callback);
            }
        }
    }

public:
    explicit CircularSupportGeneratorV2(size_t dimension) noexcept
        : dimension_(dimension), dimension_mask_(bs64::set_all_n_bits(dimension)) {}

    // Callback signature: void(bitset64 support, size_t support_size).
    template<class Callback>
    inline void generate(Callback&& callback) {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            word_.fill(0);
            emitted_ = false;
            generate_necklaces(1, 1, 0, 0, 0, callback);
            if (!emitted_)
                return;
        }
    }

    inline size_t add_forbidden(bitset64 support) {
        pending_forbidden_.push_back(support);
        // add_forbidden() is called synchronously from the callback above, so this is the multiplier calculated for that
        // emitted support.
        return emitted_multiplier_;
    }
};
