#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include <fracessa/bitset64.hpp>
#include <fracessa/bitset_multiword.hpp>

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
 * Multiword counterpart of NonCircularSupportGenerator. The search order and pruning rule are identical, but recursion mutates one
 * pre-sized mask instead of copying a vector at every branch. The callback reference is valid only until the callback returns.
 */
class NonCircularSupportGeneratorMultiword {
private:
    size_t dimension_;
    std::vector<std::vector<bitset_multiword>> forbidden_by_lowest_;
    std::vector<bitset_multiword> pending_forbidden_;
    bitset_multiword support_;
    size_t target_cardinality_ = 0;
    bool emitted_ = false;

    inline void activate_pending()
    {
        for (bitset_multiword& support : pending_forbidden_)
            forbidden_by_lowest_[support.find_pos_first_set_bit()].push_back(std::move(support));
        pending_forbidden_.clear();
    }

    inline bool completes_forbidden(size_t new_lowest_bit) const noexcept
    {
        for (const bitset_multiword& forbidden : forbidden_by_lowest_[new_lowest_bit])
            if (forbidden.is_subset_of(support_)) return true;
        return false;
    }

    template<class Callback>
    inline void generate_from(size_t bits_remaining, size_t needed, Callback& callback)
    {
        if (needed == 0) {
            emitted_ = true;
            callback(support_, target_cardinality_);
            return;
        }
        if (needed > bits_remaining) return;

        const size_t bit = bits_remaining - 1;
        if (needed < bits_remaining) generate_from(bit, needed, callback);

        support_.set_bit_at_pos(bit);
        if (!completes_forbidden(bit)) generate_from(bit, needed - 1, callback);
        support_.clear_bit_at_pos(bit);
    }

public:
    explicit NonCircularSupportGeneratorMultiword(size_t dimension)
        : dimension_(dimension), forbidden_by_lowest_(dimension), support_(dimension)
    {}

    // Callback signature: void(const bitset_multiword& support, size_t support_size).
    template<class Callback>
    inline void generate(Callback&& callback)
    {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            emitted_ = false;
            generate_from(dimension_, target_cardinality_, callback);
            if (!emitted_) return;
        }
    }

    inline void add_forbidden(const bitset_multiword& support)
    {
        assert(support.dimension() == dimension_);
        pending_forbidden_.push_back(support);
    }
};
