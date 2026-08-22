#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <fracessa/support.hpp>

namespace fracessa::support {

/** Direct fixed-density binary bracelet generator with exact forbidden-orbit pruning. */
class CircularSupportGenerator {
private:
    SupportContext& context_;
    size_t dimension_;
    size_t target_cardinality_ = 0;
    Support support_;
    Support reflected_;
    Support orbit_current_;

    uint64_t small_support_ = 0;
    std::array<size_t, kMaxBitsetDimension + 2> small_positions_{};
    std::array<uint8_t, kMaxBitsetDimension + 2> small_prefix_density_{};
    std::array<std::vector<uint64_t>, kMaxBitsetDimension> small_forbidden_by_lowest_;
    std::vector<uint64_t> small_pending_forbidden_;

    std::vector<size_t> large_positions_;
    std::vector<size_t> large_prefix_density_;
    std::vector<std::vector<Support>> large_forbidden_by_lowest_;
    std::vector<Support> large_pending_forbidden_;

    bool has_active_forbidden_ = false;
    bool emitted_ = false;

    inline void activate_pending_small()
    {
        if (!small_pending_forbidden_.empty()) has_active_forbidden_ = true;
        for (const uint64_t support : small_pending_forbidden_)
            small_forbidden_by_lowest_[ctz64(support)].push_back(support);
        small_pending_forbidden_.clear();
    }

    inline bool completes_forbidden_small(size_t new_lowest_bit) const noexcept
    {
        for (const uint64_t forbidden : small_forbidden_by_lowest_[new_lowest_bit])
            if ((forbidden & ~small_support_) == 0) return true;
        return false;
    }

    inline bool position_is_set_small(size_t position) const noexcept
    {
        return (small_support_ & (uint64_t{1} << (dimension_ - position))) != 0;
    }

    inline void set_position_small(size_t position) noexcept
    {
        small_support_ |= uint64_t{1} << (dimension_ - position);
    }

    inline void clear_position_small(size_t position) noexcept
    {
        small_support_ &= ~(uint64_t{1} << (dimension_ - position));
    }

    inline int check_reverse_small(size_t end_position) const noexcept
    {
        for (size_t position = small_positions_[1]; position <= (end_position + 1) / 2; ++position) {
            const bool forward = position_is_set_small(position);
            const bool reverse = position_is_set_small(end_position - position + 1);
            if (forward < reverse) return 1;
            if (forward > reverse) return -1;
        }
        return 0;
    }

    inline bool update_reverse_result_small(size_t density, size_t palindrome_prefix, bool reverse_smaller) const noexcept
    {
        const size_t latest_position = small_positions_[density];
        if (latest_position > (dimension_ - palindrome_prefix) / 2 + palindrome_prefix) {
            const size_t mirrored_position = dimension_ - latest_position + palindrome_prefix + 1;
            const size_t mirrored_density = small_prefix_density_[mirrored_position];
            if (mirrored_density == 0) {
                reverse_smaller = false;
            } else if (mirrored_density < density) {
                const size_t latest_zero_run = latest_position - small_positions_[density - 1] - 1;
                const size_t mirrored_zero_run = small_positions_[mirrored_density + 1] - small_positions_[mirrored_density] - 1;
                if (latest_zero_run > mirrored_zero_run) reverse_smaller = true;
            }
        }
        return reverse_smaller;
    }

    template<class Callback>
    inline void emit_final_small(size_t period_density, size_t palindrome_prefix, bool reverse_smaller, Callback& callback)
    {
        const size_t next_position = (target_cardinality_ / period_density) * small_positions_[period_density]
                                     + small_positions_[target_cardinality_ % period_density];
        if (next_position < dimension_ || (next_position == dimension_ && target_cardinality_ % period_density != 0)) return;

        set_position_small(dimension_);
        if (!completes_forbidden_small(0)) {
            reverse_smaller = update_reverse_result_small(target_cardinality_, palindrome_prefix, reverse_smaller);
            if (!reverse_smaller) {
                emitted_ = true;
                context_.set_small_bits(support_, small_support_);
                callback(support_, target_cardinality_);
            }
        }
        clear_position_small(dimension_);
    }

    template<class Callback>
    inline void generate_bracelets_small(size_t density, size_t period_density, size_t palindrome_prefix,
                                         bool reverse_smaller, Callback& callback)
    {
        reverse_smaller = update_reverse_result_small(density, palindrome_prefix, reverse_smaller);

        if (density >= target_cardinality_ - 1) {
            emit_final_small(period_density, palindrome_prefix, reverse_smaller, callback);
            return;
        }

        size_t tail = dimension_ - (target_cardinality_ - density) + 1;
        const size_t periodic_position = small_positions_[density + 1 - period_density] + small_positions_[period_density];
        if (periodic_position <= tail) {
            size_t next_palindrome_prefix = palindrome_prefix;
            bool next_reverse_smaller = reverse_smaller;
            small_positions_[density + 1] = periodic_position;
            set_position_small(periodic_position);
            small_prefix_density_[periodic_position] = static_cast<uint8_t>(density + 1);

            const size_t bit = dimension_ - periodic_position;
            if (!completes_forbidden_small(bit)) {
                bool recurse = true;
                if (small_positions_[1] == periodic_position - small_positions_[density]) {
                    const int reverse_order = check_reverse_small(periodic_position - 1);
                    if (reverse_order == 0) {
                        next_palindrome_prefix = periodic_position - 1;
                        next_reverse_smaller = false;
                    }
                    recurse = reverse_order != -1;
                }
                if (recurse)
                    generate_bracelets_small(density + 1, period_density, next_palindrome_prefix, next_reverse_smaller, callback);
            }

            clear_position_small(periodic_position);
            small_prefix_density_[periodic_position] = 0;
            tail = periodic_position - 1;
        }

        for (size_t position = tail; position > small_positions_[density]; --position) {
            small_positions_[density + 1] = position;
            set_position_small(position);
            small_prefix_density_[position] = static_cast<uint8_t>(density + 1);
            const size_t bit = dimension_ - position;
            if (!completes_forbidden_small(bit))
                generate_bracelets_small(density + 1, density + 1, palindrome_prefix, reverse_smaller, callback);
            clear_position_small(position);
            small_prefix_density_[position] = 0;
        }
    }

    template<class Callback>
    inline void generate_small(Callback& callback)
    {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending_small();
            small_support_ = 0;
            small_positions_.fill(0);
            small_prefix_density_.fill(0);
            emitted_ = false;

            if (target_cardinality_ == 1) {
                small_support_ = 1;
                if (!completes_forbidden_small(0)) {
                    emitted_ = true;
                    context_.set_small_bits(support_, small_support_);
                    callback(support_, target_cardinality_);
                }
            } else if (target_cardinality_ == dimension_) {
                if (!has_active_forbidden_) {
                    small_support_ = set_all_n_bits(dimension_);
                    emitted_ = true;
                    context_.set_small_bits(support_, small_support_);
                    callback(support_, target_cardinality_);
                }
            } else {
                small_positions_[target_cardinality_] = dimension_;
                small_prefix_density_[dimension_] = static_cast<uint8_t>(target_cardinality_);
                const size_t first_latest = dimension_ - target_cardinality_ + 1;
                const size_t first_earliest = (dimension_ - 1) / target_cardinality_ + 1;

                for (size_t first = first_latest;; --first) {
                    small_positions_[1] = first;
                    set_position_small(first);
                    small_prefix_density_[first] = 1;
                    const size_t bit = dimension_ - first;
                    if (!completes_forbidden_small(bit)) generate_bracelets_small(1, 1, first - 1, false, callback);
                    clear_position_small(first);
                    small_prefix_density_[first] = 0;
                    if (first == first_earliest) break;
                }
            }

            if (!emitted_) return;
        }
    }

    inline void activate_pending_large()
    {
        if (!large_pending_forbidden_.empty()) has_active_forbidden_ = true;
        for (Support& support : large_pending_forbidden_)
            large_forbidden_by_lowest_[context_.first(support)].push_back(std::move(support));
        large_pending_forbidden_.clear();
    }

    inline bool completes_forbidden_large(size_t new_lowest_bit) const noexcept
    {
        for (const Support& forbidden : large_forbidden_by_lowest_[new_lowest_bit])
            if (context_.is_subset_of(forbidden, support_)) return true;
        return false;
    }

    inline bool position_is_set_large(size_t position) const noexcept
    {
        return context_.contains(support_, dimension_ - position);
    }

    inline void set_position_large(size_t position) noexcept
    {
        context_.set(support_, dimension_ - position);
    }

    inline void clear_position_large(size_t position) noexcept
    {
        context_.reset(support_, dimension_ - position);
    }

    inline int check_reverse_large(size_t end_position) const noexcept
    {
        for (size_t position = large_positions_[1]; position <= (end_position + 1) / 2; ++position) {
            const bool forward = position_is_set_large(position);
            const bool reverse = position_is_set_large(end_position - position + 1);
            if (forward < reverse) return 1;
            if (forward > reverse) return -1;
        }
        return 0;
    }

    inline bool update_reverse_result_large(size_t density, size_t palindrome_prefix, bool reverse_smaller) const noexcept
    {
        const size_t latest_position = large_positions_[density];
        if (latest_position > (dimension_ - palindrome_prefix) / 2 + palindrome_prefix) {
            const size_t mirrored_position = dimension_ - latest_position + palindrome_prefix + 1;
            const size_t mirrored_density = large_prefix_density_[mirrored_position];
            if (mirrored_density == 0) {
                reverse_smaller = false;
            } else if (mirrored_density < density) {
                const size_t latest_zero_run = latest_position - large_positions_[density - 1] - 1;
                const size_t mirrored_zero_run = large_positions_[mirrored_density + 1] - large_positions_[mirrored_density] - 1;
                if (latest_zero_run > mirrored_zero_run) reverse_smaller = true;
            }
        }
        return reverse_smaller;
    }

    template<class Callback>
    inline void emit_final_large(size_t period_density, size_t palindrome_prefix, bool reverse_smaller, Callback& callback)
    {
        const size_t next_position = (target_cardinality_ / period_density) * large_positions_[period_density]
                                     + large_positions_[target_cardinality_ % period_density];
        if (next_position < dimension_ || (next_position == dimension_ && target_cardinality_ % period_density != 0)) return;

        set_position_large(dimension_);
        if (!completes_forbidden_large(0)) {
            reverse_smaller = update_reverse_result_large(target_cardinality_, palindrome_prefix, reverse_smaller);
            if (!reverse_smaller) {
                emitted_ = true;
                callback(support_, target_cardinality_);
            }
        }
        clear_position_large(dimension_);
    }

    template<class Callback>
    inline void generate_bracelets_large(size_t density, size_t period_density, size_t palindrome_prefix,
                                         bool reverse_smaller, Callback& callback)
    {
        reverse_smaller = update_reverse_result_large(density, palindrome_prefix, reverse_smaller);

        if (density >= target_cardinality_ - 1) {
            emit_final_large(period_density, palindrome_prefix, reverse_smaller, callback);
            return;
        }

        size_t tail = dimension_ - (target_cardinality_ - density) + 1;
        const size_t periodic_position = large_positions_[density + 1 - period_density] + large_positions_[period_density];
        if (periodic_position <= tail) {
            size_t next_palindrome_prefix = palindrome_prefix;
            bool next_reverse_smaller = reverse_smaller;
            large_positions_[density + 1] = periodic_position;
            set_position_large(periodic_position);
            large_prefix_density_[periodic_position] = density + 1;

            const size_t bit = dimension_ - periodic_position;
            if (!completes_forbidden_large(bit)) {
                bool recurse = true;
                if (large_positions_[1] == periodic_position - large_positions_[density]) {
                    const int reverse_order = check_reverse_large(periodic_position - 1);
                    if (reverse_order == 0) {
                        next_palindrome_prefix = periodic_position - 1;
                        next_reverse_smaller = false;
                    }
                    recurse = reverse_order != -1;
                }
                if (recurse)
                    generate_bracelets_large(density + 1, period_density, next_palindrome_prefix, next_reverse_smaller, callback);
            }

            clear_position_large(periodic_position);
            large_prefix_density_[periodic_position] = 0;
            tail = periodic_position - 1;
        }

        for (size_t position = tail; position > large_positions_[density]; --position) {
            large_positions_[density + 1] = position;
            set_position_large(position);
            large_prefix_density_[position] = density + 1;
            const size_t bit = dimension_ - position;
            if (!completes_forbidden_large(bit))
                generate_bracelets_large(density + 1, density + 1, palindrome_prefix, reverse_smaller, callback);
            clear_position_large(position);
            large_prefix_density_[position] = 0;
        }
    }

    template<class Callback>
    inline void generate_large(Callback& callback)
    {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending_large();
            context_.clear(support_);
            std::fill(large_positions_.begin(), large_positions_.end(), size_t{0});
            std::fill(large_prefix_density_.begin(), large_prefix_density_.end(), size_t{0});
            emitted_ = false;

            if (target_cardinality_ == 1) {
                set_position_large(dimension_);
                if (!completes_forbidden_large(0)) {
                    emitted_ = true;
                    callback(support_, target_cardinality_);
                }
                clear_position_large(dimension_);
            } else if (target_cardinality_ == dimension_) {
                if (!has_active_forbidden_) {
                    context_.set_all(support_);
                    emitted_ = true;
                    callback(support_, target_cardinality_);
                    context_.clear(support_);
                }
            } else {
                large_positions_[target_cardinality_] = dimension_;
                large_prefix_density_[dimension_] = target_cardinality_;
                const size_t first_latest = dimension_ - target_cardinality_ + 1;
                const size_t first_earliest = (dimension_ - 1) / target_cardinality_ + 1;

                for (size_t first = first_latest;; --first) {
                    large_positions_[1] = first;
                    set_position_large(first);
                    large_prefix_density_[first] = 1;
                    const size_t bit = dimension_ - first;
                    if (!completes_forbidden_large(bit)) generate_bracelets_large(1, 1, first - 1, false, callback);
                    clear_position_large(first);
                    large_prefix_density_[first] = 0;
                    if (first == first_earliest) break;
                }
            }

            if (!emitted_) return;
        }
    }

public:
    explicit CircularSupportGenerator(SupportContext& context)
        : context_(context)
        , dimension_(context.dimension())
        , support_(context.make())
        , reflected_(context.make())
        , orbit_current_(context.make())
    {
        if (!context_.is_small()) {
            large_positions_.assign(dimension_ + 2, 0);
            large_prefix_density_.assign(dimension_ + 2, 0);
            large_forbidden_by_lowest_.resize(dimension_);
        }
    }

    template<class Callback>
    inline void generate(Callback&& callback)
    {
        if (context_.is_small()) generate_small(callback);
        else generate_large(callback);
    }

    inline size_t add_forbidden_orbit(const Support& support)
    {
        if (context_.is_small()) {
            const uint64_t support_bits = context_.small_bits(support);
            const uint64_t reflected = reflect(support_bits, dimension_);
            bool reflection_is_rotation = false;
            size_t multiplier = 0;
            uint64_t current = support_bits;
            do {
                reflection_is_rotation = reflection_is_rotation || current == reflected;
                small_pending_forbidden_.push_back(current);
                ++multiplier;
                current = rot_one_right(current, dimension_);
            } while (current != support_bits);

            if (!reflection_is_rotation) {
                current = reflected;
                do {
                    small_pending_forbidden_.push_back(current);
                    ++multiplier;
                    current = rot_one_right(current, dimension_);
                } while (current != reflected);
            }
            return multiplier;
        }

        context_.copy(reflected_, support);
        context_.reflect(reflected_);

        bool reflection_is_rotation = false;
        size_t multiplier = 0;
        context_.copy(orbit_current_, support);
        do {
            reflection_is_rotation = reflection_is_rotation || context_.equal(orbit_current_, reflected_);
            large_pending_forbidden_.push_back(context_.clone(orbit_current_));
            ++multiplier;
            context_.rotate_one_right(orbit_current_);
        } while (!context_.equal(orbit_current_, support));

        if (!reflection_is_rotation) {
            context_.copy(orbit_current_, reflected_);
            do {
                large_pending_forbidden_.push_back(context_.clone(orbit_current_));
                ++multiplier;
                context_.rotate_one_right(orbit_current_);
            } while (!context_.equal(orbit_current_, reflected_));
        }
        return multiplier;
    }

    inline bool has_supports_after_singletons()
    {
        assert(target_cardinality_ == 1);
        if (context_.is_small()) {
            uint64_t forbidden_singletons = 0;
            for (const uint64_t support : small_pending_forbidden_) {
                const size_t cardinality = popcount64(support);
                assert(cardinality == 1);
                if (cardinality == 1) forbidden_singletons |= support;
            }
            return popcount64(forbidden_singletons) + 1 < dimension_;
        }

        Support forbidden_singletons = context_.make();
        for (const Support& support : large_pending_forbidden_) {
            const size_t cardinality = context_.count(support);
            assert(cardinality == 1);
            if (cardinality == 1) context_.add(forbidden_singletons, support);
        }
        return context_.count(forbidden_singletons) + 1 < dimension_;
    }
};

} // namespace fracessa::support
