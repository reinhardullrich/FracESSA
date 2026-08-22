#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include <fracessa/bitset.hpp>
#include <fracessa/bitset_multiword.hpp>

namespace fracessa::support {

/**
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
    struct invader_rule {
        bitset required;
        bitset invaders;
    };

    size_t dimension_;
    // A rule is checked when its lowest set bit is added, because that is the first moment the partial support can
    // contain the complete rule.
    std::array<std::vector<bitset>, kMaxBitsetDimension> forbidden_by_lowest_;
    // Rules found by the callback stay pending until the next support size.
    std::vector<bitset> pending_forbidden_;
    std::array<std::vector<invader_rule>, kMaxBitsetDimension> invader_by_lowest_;
    std::vector<invader_rule> pending_invader_;
    std::vector<const invader_rule*> active_invader_;
    size_t target_cardinality_ = 0;
    bool emitted_ = false;

    inline void activate_pending() {
        for (const bitset support : pending_forbidden_)
            forbidden_by_lowest_[ctz64(support)].push_back(support);
        pending_forbidden_.clear();
        for (const invader_rule& rule : pending_invader_)
            invader_by_lowest_[ctz64(rule.required)].push_back(rule);
        pending_invader_.clear();
    }

    inline bool activate_invader_rules(bitset partial, size_t new_lowest_bit) {
        const bitset future = new_lowest_bit == 0 ? 0 : (bitset{1} << new_lowest_bit) - 1;
        for (const invader_rule& rule : invader_by_lowest_[new_lowest_bit]) {
            if (!is_subset_of(rule.required, partial) || (rule.invaders & partial) != 0) continue;
            if ((rule.invaders & future) == 0) return true;
            active_invader_.push_back(&rule);
        }
        return false;
    }

    inline bool violates_invader_rule(bitset support) const noexcept {
        for (const invader_rule* rule : active_invader_)
            if ((rule->invaders & support) == 0) return true;
        return false;
    }

    inline bool completes_forbidden(bitset partial, size_t new_lowest_bit) const noexcept {
        for (const bitset forbidden : forbidden_by_lowest_[new_lowest_bit]) {
            if (is_subset_of(forbidden, partial))
                return true;
        }
        return false;
    }

    template<class Callback>
    inline void generate_from(size_t bits_remaining, size_t needed, bitset partial, Callback& callback) {
        if (needed == 0) {
            if (violates_invader_rule(partial)) return;
            emitted_ = true;
            callback(partial, target_cardinality_);
            return;
        }
        if (needed > bits_remaining)
            return;

        const size_t bit = bits_remaining - 1;
        if (needed < bits_remaining)
            generate_from(bit, needed, partial, callback);

        const bitset with_bit = set_bit_at_pos(partial, bit);
        if (!completes_forbidden(with_bit, bit)) {
            const size_t active_size = active_invader_.size();
            if (!activate_invader_rules(with_bit, bit)) generate_from(bit, needed - 1, with_bit, callback);
            active_invader_.resize(active_size);
        }
    }

public:
    explicit NonCircularSupportGenerator(size_t dimension) noexcept : dimension_(dimension) {}

    // Callback signature: void(bitset support, size_t support_size).
    // It is called synchronously once for each generated support.
    template<class Callback>
    inline void generate(Callback&& callback) {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            assert(active_invader_.empty());
            emitted_ = false;
            generate_from(dimension_, target_cardinality_, 0, callback);
            // If no support of this size survives, no larger support can survive either: every larger support contains
            // one of this size.
            if (!emitted_)
                return;
        }
    }

    inline void add_forbidden(bitset support) {
        pending_forbidden_.push_back(support);
    }

    inline void add_invader_clause(bitset required, bitset invaders) {
        assert(required != 0 && invaders != 0 && (required & invaders) == 0);
        pending_invader_.push_back({required, invaders});
    }

    // Call from the final singleton callback. Pending rules then contain exactly the singleton roots found in that layer.
    inline bool has_supports_after_singletons() const noexcept {
        assert(target_cardinality_ == 1);
        bitset forbidden_singletons = 0;
        for (const bitset support : pending_forbidden_) {
            const size_t cardinality = count_set_bits(support);
            assert(cardinality == 1);
            if (cardinality == 1) forbidden_singletons |= support;
        }
        return count_set_bits(forbidden_singletons) + 1 < dimension_;
    }
};

/**
 * Multiword counterpart of NonCircularSupportGenerator. The search order and pruning rule are identical, but recursion mutates one
 * pre-sized mask instead of copying a vector at every branch. The callback reference is valid only until the callback returns.
 */
class NonCircularSupportGeneratorMultiword {
private:
    struct invader_rule {
        bitset_multiword required;
        std::vector<size_t> invaders;

        invader_rule(const bitset_multiword& required_mask, const bitset_multiword& invader_mask)
            : required(required_mask)
        {
            invader_mask.extract_set_indices(invaders);
        }
    };

    size_t dimension_;
    std::vector<std::vector<bitset_multiword>> forbidden_by_lowest_;
    std::vector<bitset_multiword> pending_forbidden_;
    bitset_multiword support_;
    std::vector<std::vector<invader_rule>> invader_by_lowest_;
    std::vector<invader_rule> pending_invader_;
    std::vector<const invader_rule*> active_invader_;
    size_t target_cardinality_ = 0;
    bool emitted_ = false;

    inline void activate_pending()
    {
        for (bitset_multiword& support : pending_forbidden_)
            forbidden_by_lowest_[support.find_pos_first_set_bit()].push_back(std::move(support));
        pending_forbidden_.clear();
        for (invader_rule& rule : pending_invader_)
            invader_by_lowest_[rule.required.find_pos_first_set_bit()].push_back(std::move(rule));
        pending_invader_.clear();
    }

    inline bool invader_selected(const invader_rule& rule) const noexcept
    {
        for (const size_t position : rule.invaders)
            if (support_.is_set_at_pos(position)) return true;
        return false;
    }

    inline bool activate_invader_rules(size_t new_lowest_bit)
    {
        for (const invader_rule& rule : invader_by_lowest_[new_lowest_bit]) {
            if (!rule.required.is_subset_of(support_) || invader_selected(rule)) continue;
            bool has_future_invader = false;
            for (const size_t position : rule.invaders) has_future_invader |= position < new_lowest_bit;
            if (!has_future_invader) return true;
            active_invader_.push_back(&rule);
        }
        return false;
    }

    inline bool violates_invader_rule() const noexcept
    {
        for (const invader_rule* rule : active_invader_)
            if (!invader_selected(*rule)) return true;
        return false;
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
            if (violates_invader_rule()) return;
            emitted_ = true;
            callback(support_, target_cardinality_);
            return;
        }
        if (needed > bits_remaining) return;

        const size_t bit = bits_remaining - 1;
        if (needed < bits_remaining) generate_from(bit, needed, callback);

        support_.set_bit_at_pos(bit);
        if (!completes_forbidden(bit)) {
            const size_t active_size = active_invader_.size();
            if (!activate_invader_rules(bit)) generate_from(bit, needed - 1, callback);
            active_invader_.resize(active_size);
        }
        support_.clear_bit_at_pos(bit);
    }

public:
    explicit NonCircularSupportGeneratorMultiword(size_t dimension)
        : dimension_(dimension), forbidden_by_lowest_(dimension), support_(dimension), invader_by_lowest_(dimension)
    {}

    // Callback signature: void(const bitset_multiword& support, size_t support_size).
    template<class Callback>
    inline void generate(Callback&& callback)
    {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            assert(active_invader_.empty());
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

    inline void add_invader_clause(const bitset_multiword& required, const bitset_multiword& invaders)
    {
        assert(required.dimension() == dimension_ && invaders.dimension() == dimension_ && !invaders.empty());
        pending_invader_.emplace_back(required, invaders);
    }

    // Call from the final singleton callback. Pending rules then contain exactly the singleton roots found in that layer.
    inline bool has_supports_after_singletons() const
    {
        assert(target_cardinality_ == 1);
        bitset_multiword forbidden_singletons(dimension_);
        for (const bitset_multiword& support : pending_forbidden_) {
            const size_t cardinality = support.count_set_bits();
            assert(cardinality == 1);
            if (cardinality == 1) forbidden_singletons.set_bit_at_pos(support.find_pos_first_set_bit());
        }
        return forbidden_singletons.count_set_bits() + 1 < dimension_;
    }
};

} // namespace fracessa::support
