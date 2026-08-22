#pragma once

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include <fracessa/support.hpp>

namespace fracessa::support {

/** Generate every nonempty non-circular support in increasing cardinality and numeric-mask order. */
class NonCircularSupportGenerator {
private:
    SupportContext& context_;
    std::vector<std::vector<Support>> forbidden_by_lowest_;
    std::vector<Support> pending_forbidden_;
    Support support_;
    size_t target_cardinality_ = 0;
    bool emitted_ = false;

    inline void activate_pending()
    {
        for (Support& support : pending_forbidden_)
            forbidden_by_lowest_[context_.first(support)].push_back(std::move(support));
        pending_forbidden_.clear();
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
            emitted_ = true;
            callback(support_, target_cardinality_);
            return;
        }
        if (needed > bits_remaining) return;

        const size_t bit = bits_remaining - 1;
        if (needed < bits_remaining) generate_from(bit, needed, callback);

        context_.set(support_, bit);
        if (!completes_forbidden(bit)) generate_from(bit, needed - 1, callback);
        context_.reset(support_, bit);
    }

public:
    explicit NonCircularSupportGenerator(SupportContext& context)
        : context_(context)
        , forbidden_by_lowest_(context.dimension())
        , support_(context.make())
    {}

    // The callback reference is valid only until the callback returns.
    template<class Callback>
    inline void generate(Callback&& callback)
    {
        for (target_cardinality_ = 1; target_cardinality_ <= context_.dimension(); ++target_cardinality_) {
            activate_pending();
            emitted_ = false;
            generate_from(context_.dimension(), target_cardinality_, callback);
            if (!emitted_) return;
        }
    }

    inline void add_forbidden(const Support& support)
    {
        pending_forbidden_.push_back(context_.clone(support));
    }

    // Call from the final singleton callback.
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
