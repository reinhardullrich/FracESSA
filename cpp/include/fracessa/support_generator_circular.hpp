#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <fracessa/bitset.hpp>
#include <fracessa/bitset_multiword.hpp>

namespace fracessa::support {

/*
 * Direct fixed-density bracelet generation from Karim, Alamgir, and Husnine (2014), specialized to binary supports.
 * A bracelet is a class of supports that differ only by rotation or reflection. For each support size, this generator emits the
 * numerically smallest support in each class. It compares the growing word with its reflection during recursion, so a mirrored copy
 * is rejected before it is fully constructed. The independent reference generator in the test tree verifies the callback, pruning,
 * and orbit-multiplier contract.
 */
class CircularSupportGenerator {
private:
    size_t dimension_;
    size_t target_cardinality_ = 0;
    bitset support_ = 0;
    // positions_[density] is the one-indexed position of that density's 1-bit. prefix_density_ maps selected positions
    // back to their density, which makes the paper's reversal comparison constant-time between successive 1-bits.
    std::array<size_t, kMaxBitsetDimension + 2> positions_{};
    std::array<uint8_t, kMaxBitsetDimension + 2> prefix_density_{};
    std::array<std::vector<bitset>, kMaxBitsetDimension> forbidden_by_lowest_;
    std::vector<bitset> pending_forbidden_;
    bool has_active_forbidden_ = false;
    bool emitted_ = false;

    inline void activate_pending() {
        if (!pending_forbidden_.empty())
            has_active_forbidden_ = true;
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

    inline bool position_is_set(size_t position) const noexcept {
        return (support_ & (1ULL << (dimension_ - position))) != 0;
    }

    inline void set_position(size_t position) noexcept {
        support_ |= 1ULL << (dimension_ - position);
    }

    inline void clear_position(size_t position) noexcept {
        support_ &= ~(1ULL << (dimension_ - position));
    }

    // Return 1 when the current prefix is smaller than its reversal, 0 when equal, and -1 when the reversal is smaller.
    inline int check_reverse(size_t end_position) const noexcept {
        for (size_t position = positions_[1]; position <= (end_position + 1) / 2; ++position) {
            const bool forward = position_is_set(position);
            const bool reverse = position_is_set(end_position - position + 1);
            if (forward < reverse)
                return 1;
            if (forward > reverse)
                return -1;
        }
        return 0;
    }

    inline bool update_reverse_result(size_t density, size_t palindrome_prefix, bool reverse_smaller) const noexcept {
        const size_t latest_position = positions_[density];
        if (latest_position > (dimension_ - palindrome_prefix) / 2 + palindrome_prefix) {
            const size_t mirrored_position = dimension_ - latest_position + palindrome_prefix + 1;
            const size_t mirrored_density = prefix_density_[mirrored_position];
            // Every newly placed binary symbol is 1, so a zero mirrored density means 1 > 0.
            if (mirrored_density == 0) {
                reverse_smaller = false;
            } else if (mirrored_density < density) {
                // Equality identifies the latest symbol itself, whose following zero run is not known yet.
                const size_t latest_zero_run = latest_position - positions_[density - 1] - 1;
                const size_t mirrored_zero_run = positions_[mirrored_density + 1] - positions_[mirrored_density] - 1;
                if (latest_zero_run > mirrored_zero_run)
                    reverse_smaller = true;
            }
        }
        return reverse_smaller;
    }

    template<class Callback>
    inline void emit_final(size_t period_density, size_t palindrome_prefix, bool reverse_smaller, Callback& callback) {
        // In the paper's PrintD step, these are exactly the two cases where the final 1-bit belongs at position n.
        const size_t next_position = (target_cardinality_ / period_density) * positions_[period_density]
                                     + positions_[target_cardinality_ % period_density];
        if (next_position < dimension_
            || (next_position == dimension_ && target_cardinality_ % period_density != 0))
            return;

        set_position(dimension_);
        if (!completes_forbidden(support_, 0)) {
            reverse_smaller = update_reverse_result(target_cardinality_, palindrome_prefix, reverse_smaller);
            if (!reverse_smaller) {
                emitted_ = true;
                callback(support_, target_cardinality_);
            }
        }
        clear_position(dimension_);
    }

    template<class Callback>
    inline void generate_bracelets(size_t density, size_t period_density, size_t palindrome_prefix,
                                   bool reverse_smaller, Callback& callback) {
        reverse_smaller = update_reverse_result(density, palindrome_prefix, reverse_smaller);

        if (density >= target_cardinality_ - 1) {
            emit_final(period_density, palindrome_prefix, reverse_smaller, callback);
            return;
        }

        size_t tail = dimension_ - (target_cardinality_ - density) + 1;
        const size_t periodic_position = positions_[density + 1 - period_density] + positions_[period_density];
        if (periodic_position <= tail) {
            size_t next_palindrome_prefix = palindrome_prefix;
            bool next_reverse_smaller = reverse_smaller;
            positions_[density + 1] = periodic_position;
            set_position(periodic_position);
            prefix_density_[periodic_position] = static_cast<uint8_t>(density + 1);

            const size_t bit = dimension_ - periodic_position;
            if (!completes_forbidden(support_, bit)) {
                bool recurse = true;
                if (positions_[1] == periodic_position - positions_[density]) {
                    const int reverse_order = check_reverse(periodic_position - 1);
                    if (reverse_order == 0) {
                        next_palindrome_prefix = periodic_position - 1;
                        next_reverse_smaller = false;
                    }
                    recurse = reverse_order != -1;
                }
                if (recurse)
                    generate_bracelets(density + 1, period_density, next_palindrome_prefix, next_reverse_smaller, callback);
            }

            clear_position(periodic_position);
            prefix_density_[periodic_position] = 0;
            tail = periodic_position - 1;
        }

        for (size_t position = tail; position > positions_[density]; --position) {
            positions_[density + 1] = position;
            set_position(position);
            prefix_density_[position] = static_cast<uint8_t>(density + 1);
            const size_t bit = dimension_ - position;
            if (!completes_forbidden(support_, bit))
                generate_bracelets(density + 1, density + 1, palindrome_prefix, reverse_smaller, callback);
            clear_position(position);
            prefix_density_[position] = 0;
        }
    }

public:
    explicit CircularSupportGenerator(size_t dimension) noexcept : dimension_(dimension) {}

    // Callback signature: void(bitset support, size_t support_size).
    template<class Callback>
    inline void generate(Callback&& callback) {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            support_ = 0;
            positions_.fill(0);
            prefix_density_.fill(0);
            emitted_ = false;

            if (target_cardinality_ == 1) {
                support_ = 1;
                if (!completes_forbidden(support_, 0)) {
                    emitted_ = true;
                    callback(support_, target_cardinality_);
                }
            } else if (target_cardinality_ == dimension_) {
                // The full support contains every nonempty active forbidden support.
                if (!has_active_forbidden_) {
                    support_ = set_all_n_bits(dimension_);
                    emitted_ = true;
                    callback(support_, target_cardinality_);
                }
            } else {
                positions_[target_cardinality_] = dimension_;
                prefix_density_[dimension_] = static_cast<uint8_t>(target_cardinality_);
                const size_t first_latest = dimension_ - target_cardinality_ + 1;
                const size_t first_earliest = (dimension_ - 1) / target_cardinality_ + 1;

                for (size_t first = first_latest;; --first) {
                    positions_[1] = first;
                    set_position(first);
                    prefix_density_[first] = 1;
                    const size_t bit = dimension_ - first;
                    if (!completes_forbidden(support_, bit))
                        generate_bracelets(1, 1, first - 1, false, callback);
                    clear_position(first);
                    prefix_density_[first] = 0;
                    if (first == first_earliest)
                        break;
                }
            }

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

/*
 * Multiword counterpart of CircularSupportGenerator. It retains the direct fixed-density bracelet algorithm and mutates one
 * pre-sized support mask during recursion. The callback reference is valid only until the callback returns.
 */
class CircularSupportGeneratorMultiword {
private:
    size_t dimension_;
    size_t target_cardinality_ = 0;
    bitset_multiword support_;
    bitset_multiword reflected_;
    bitset_multiword orbit_current_;
    std::vector<size_t> positions_;
    std::vector<size_t> prefix_density_;
    std::vector<std::vector<bitset_multiword>> forbidden_by_lowest_;
    std::vector<bitset_multiword> pending_forbidden_;
    bool has_active_forbidden_ = false;
    bool emitted_ = false;

    inline void activate_pending()
    {
        if (!pending_forbidden_.empty()) has_active_forbidden_ = true;
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

    inline bool position_is_set(size_t position) const noexcept
    {
        return support_.is_set_at_pos(dimension_ - position);
    }

    inline void set_position(size_t position) noexcept
    {
        support_.set_bit_at_pos(dimension_ - position);
    }

    inline void clear_position(size_t position) noexcept
    {
        support_.clear_bit_at_pos(dimension_ - position);
    }

    // Return 1 when the current prefix is smaller than its reversal, 0 when equal, and -1 when the reversal is smaller.
    inline int check_reverse(size_t end_position) const noexcept
    {
        for (size_t position = positions_[1]; position <= (end_position + 1) / 2; ++position) {
            const bool forward = position_is_set(position);
            const bool reverse = position_is_set(end_position - position + 1);
            if (forward < reverse) return 1;
            if (forward > reverse) return -1;
        }
        return 0;
    }

    inline bool update_reverse_result(size_t density, size_t palindrome_prefix, bool reverse_smaller) const noexcept
    {
        const size_t latest_position = positions_[density];
        if (latest_position > (dimension_ - palindrome_prefix) / 2 + palindrome_prefix) {
            const size_t mirrored_position = dimension_ - latest_position + palindrome_prefix + 1;
            const size_t mirrored_density = prefix_density_[mirrored_position];
            if (mirrored_density == 0) {
                reverse_smaller = false;
            } else if (mirrored_density < density) {
                const size_t latest_zero_run = latest_position - positions_[density - 1] - 1;
                const size_t mirrored_zero_run = positions_[mirrored_density + 1] - positions_[mirrored_density] - 1;
                if (latest_zero_run > mirrored_zero_run) reverse_smaller = true;
            }
        }
        return reverse_smaller;
    }

    template<class Callback>
    inline void emit_final(size_t period_density, size_t palindrome_prefix, bool reverse_smaller, Callback& callback)
    {
        const size_t next_position = (target_cardinality_ / period_density) * positions_[period_density]
                                     + positions_[target_cardinality_ % period_density];
        if (next_position < dimension_ || (next_position == dimension_ && target_cardinality_ % period_density != 0)) return;

        set_position(dimension_);
        if (!completes_forbidden(0)) {
            reverse_smaller = update_reverse_result(target_cardinality_, palindrome_prefix, reverse_smaller);
            if (!reverse_smaller) {
                emitted_ = true;
                callback(support_, target_cardinality_);
            }
        }
        clear_position(dimension_);
    }

    template<class Callback>
    inline void generate_bracelets(size_t density, size_t period_density, size_t palindrome_prefix,
                                   bool reverse_smaller, Callback& callback)
    {
        reverse_smaller = update_reverse_result(density, palindrome_prefix, reverse_smaller);

        if (density >= target_cardinality_ - 1) {
            emit_final(period_density, palindrome_prefix, reverse_smaller, callback);
            return;
        }

        size_t tail = dimension_ - (target_cardinality_ - density) + 1;
        const size_t periodic_position = positions_[density + 1 - period_density] + positions_[period_density];
        if (periodic_position <= tail) {
            size_t next_palindrome_prefix = palindrome_prefix;
            bool next_reverse_smaller = reverse_smaller;
            positions_[density + 1] = periodic_position;
            set_position(periodic_position);
            prefix_density_[periodic_position] = density + 1;

            const size_t bit = dimension_ - periodic_position;
            if (!completes_forbidden(bit)) {
                bool recurse = true;
                if (positions_[1] == periodic_position - positions_[density]) {
                    const int reverse_order = check_reverse(periodic_position - 1);
                    if (reverse_order == 0) {
                        next_palindrome_prefix = periodic_position - 1;
                        next_reverse_smaller = false;
                    }
                    recurse = reverse_order != -1;
                }
                if (recurse)
                    generate_bracelets(density + 1, period_density, next_palindrome_prefix, next_reverse_smaller, callback);
            }

            clear_position(periodic_position);
            prefix_density_[periodic_position] = 0;
            tail = periodic_position - 1;
        }

        for (size_t position = tail; position > positions_[density]; --position) {
            positions_[density + 1] = position;
            set_position(position);
            prefix_density_[position] = density + 1;
            const size_t bit = dimension_ - position;
            if (!completes_forbidden(bit))
                generate_bracelets(density + 1, density + 1, palindrome_prefix, reverse_smaller, callback);
            clear_position(position);
            prefix_density_[position] = 0;
        }
    }

public:
    explicit CircularSupportGeneratorMultiword(size_t dimension)
        : dimension_(dimension)
        , support_(dimension)
        , reflected_(dimension)
        , orbit_current_(dimension)
        , positions_(dimension + 2, 0)
        , prefix_density_(dimension + 2, 0)
        , forbidden_by_lowest_(dimension)
    {}

    // Callback signature: void(const bitset_multiword& support, size_t support_size).
    template<class Callback>
    inline void generate(Callback&& callback)
    {
        for (target_cardinality_ = 1; target_cardinality_ <= dimension_; ++target_cardinality_) {
            activate_pending();
            support_.clear_all();
            std::fill(positions_.begin(), positions_.end(), size_t{0});
            std::fill(prefix_density_.begin(), prefix_density_.end(), size_t{0});
            emitted_ = false;

            if (target_cardinality_ == 1) {
                set_position(dimension_);
                if (!completes_forbidden(0)) {
                    emitted_ = true;
                    callback(support_, target_cardinality_);
                }
                clear_position(dimension_);
            } else if (target_cardinality_ == dimension_) {
                if (!has_active_forbidden_) {
                    support_.set_all();
                    emitted_ = true;
                    callback(support_, target_cardinality_);
                    support_.clear_all();
                }
            } else {
                positions_[target_cardinality_] = dimension_;
                prefix_density_[dimension_] = target_cardinality_;
                const size_t first_latest = dimension_ - target_cardinality_ + 1;
                const size_t first_earliest = (dimension_ - 1) / target_cardinality_ + 1;

                for (size_t first = first_latest;; --first) {
                    positions_[1] = first;
                    set_position(first);
                    prefix_density_[first] = 1;
                    const size_t bit = dimension_ - first;
                    if (!completes_forbidden(bit)) generate_bracelets(1, 1, first - 1, false, callback);
                    clear_position(first);
                    prefix_density_[first] = 0;
                    if (first == first_earliest) break;
                }
            }

            if (!emitted_) return;
        }
    }

    inline size_t add_forbidden_orbit(const bitset_multiword& support)
    {
        assert(support.dimension() == dimension_);
        reflected_.copy_from(support);
        reflected_.reflect();

        bool reflection_is_rotation = false;
        size_t multiplier = 0;
        orbit_current_.copy_from(support);
        do {
            reflection_is_rotation = reflection_is_rotation || orbit_current_ == reflected_;
            pending_forbidden_.push_back(orbit_current_);
            ++multiplier;
            orbit_current_.rotate_one_right();
        } while (orbit_current_ != support);

        if (!reflection_is_rotation) {
            orbit_current_.copy_from(reflected_);
            do {
                pending_forbidden_.push_back(orbit_current_);
                ++multiplier;
                orbit_current_.rotate_one_right();
            } while (orbit_current_ != reflected_);
        }
        return multiplier;
    }
};

} // namespace fracessa::support
