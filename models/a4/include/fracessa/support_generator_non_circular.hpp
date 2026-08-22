#pragma once

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include <fracessa/support.hpp>

namespace fracessa::support {

/** Non-circular support generation with A4's exact sideways invader clauses. */
class NonCircularSupportGenerator {
private:
    struct invader_rule {
        Support required;
        std::vector<size_t> invaders;

        invader_rule(SupportContext& context, const Support& required_mask, const Support& invader_mask)
            : required(context.clone(required_mask))
        {
            context.extract_set_indices(invader_mask, invaders);
        }

        invader_rule(invader_rule&&) noexcept = default;
        invader_rule& operator=(invader_rule&&) noexcept = default;
        invader_rule(const invader_rule&) = delete;
        invader_rule& operator=(const invader_rule&) = delete;
    };

    SupportContext& context_;
    std::vector<std::vector<Support>> forbidden_by_lowest_;
    std::vector<Support> pending_forbidden_;
    Support support_;
    std::vector<std::vector<invader_rule>> invader_by_lowest_;
    std::vector<invader_rule> pending_invader_;
    std::vector<const invader_rule*> active_invader_;
    size_t target_cardinality_ = 0;
    bool emitted_ = false;

    inline void activate_pending()
    {
        for (Support& support : pending_forbidden_)
            forbidden_by_lowest_[context_.first(support)].push_back(std::move(support));
        pending_forbidden_.clear();
        for (invader_rule& rule : pending_invader_)
            invader_by_lowest_[context_.first(rule.required)].push_back(std::move(rule));
        pending_invader_.clear();
    }

    inline bool invader_selected(const invader_rule& rule) const noexcept
    {
        for (const size_t position : rule.invaders)
            if (context_.contains(support_, position)) return true;
        return false;
    }

    inline bool activate_invader_rules(size_t new_lowest_bit)
    {
        for (const invader_rule& rule : invader_by_lowest_[new_lowest_bit]) {
            if (!context_.is_subset_of(rule.required, support_) || invader_selected(rule)) continue;

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
        for (const Support& forbidden : forbidden_by_lowest_[new_lowest_bit])
            if (context_.is_subset_of(forbidden, support_)) return true;
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

        context_.set(support_, bit);
        if (!completes_forbidden(bit)) {
            const size_t active_size = active_invader_.size();
            if (!activate_invader_rules(bit)) generate_from(bit, needed - 1, callback);
            active_invader_.resize(active_size);
        }
        context_.reset(support_, bit);
    }

public:
    explicit NonCircularSupportGenerator(SupportContext& context)
        : context_(context)
        , forbidden_by_lowest_(context.dimension())
        , support_(context.make())
        , invader_by_lowest_(context.dimension())
    {}

    template<class Callback>
    inline void generate(Callback&& callback)
    {
        for (target_cardinality_ = 1; target_cardinality_ <= context_.dimension(); ++target_cardinality_) {
            activate_pending();
            assert(active_invader_.empty());
            emitted_ = false;
            generate_from(context_.dimension(), target_cardinality_, callback);
            if (!emitted_) return;
        }
    }

    inline void add_forbidden(const Support& support)
    {
        pending_forbidden_.push_back(context_.clone(support));
    }

    inline void add_invader_clause(const Support& required, const Support& invaders)
    {
        assert(!context_.empty(required) && !context_.empty(invaders));
        pending_invader_.emplace_back(context_, required, invaders);
    }

    inline bool has_supports_after_singletons()
    {
        assert(target_cardinality_ == 1);
        Support forbidden_singletons = context_.make();
        for (const Support& support : pending_forbidden_) {
            const size_t cardinality = context_.count(support);
            assert(cardinality == 1);
            if (cardinality == 1) context_.add(forbidden_singletons, support);
        }
        return context_.count(forbidden_singletons) + 1 < context_.dimension();
    }
};

} // namespace fracessa::support
