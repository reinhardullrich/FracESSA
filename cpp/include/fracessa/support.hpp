#pragma once

#include <fracessa/bitset.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fracessa::support {

class SupportContext;

/** One support mask. Its interpretation belongs to the matrix-wide SupportContext. */
class Support {
    friend class SupportContext;

public:
    Support(const Support&) = delete;
    Support& operator=(const Support&) = delete;

    Support(Support&& other) noexcept : storage_(std::exchange(other.storage_, 0)) {}

    Support& operator=(Support&& other) noexcept
    {
        storage_ = std::exchange(other.storage_, 0);
        return *this;
    }

private:
    explicit Support(uint64_t storage) noexcept : storage_(storage) {}

    uint64_t storage_;
};

static_assert(sizeof(Support) == sizeof(uint64_t));
static_assert(sizeof(std::uintptr_t) <= sizeof(uint64_t));

/**
 * Matrix-wide support representation and storage.
 *
 * Through dimension 64, a Support stores its bits directly. Above dimension 64, it stores a pointer to an exact-sized word array
 * owned by this context. The dimension decision and allocation bookkeeping therefore exist once per analyzer, not in every mask.
 */
class SupportContext {
public:
    explicit SupportContext(size_t dimension)
        : dimension_(dimension)
        , word_count_(dimension / 64 + static_cast<size_t>(dimension % 64 != 0))
        , last_word_mask_(make_last_word_mask(dimension))
    {
        if (dimension == 0) throw std::invalid_argument("SupportContext dimension must be positive");
    }

    SupportContext(const SupportContext&) = delete;
    SupportContext& operator=(const SupportContext&) = delete;
    SupportContext(SupportContext&&) = delete;
    SupportContext& operator=(SupportContext&&) = delete;

    size_t dimension() const noexcept { return dimension_; }
    bool is_small() const noexcept { return dimension_ <= kMaxBitsetDimension; }

    uint64_t small_bits(const Support& value) const noexcept
    {
        assert(is_small());
        return value.storage_;
    }

    void set_small_bits(Support& value, uint64_t bits) const noexcept
    {
        assert(is_small());
        value.storage_ = bits;
    }

    Support make()
    {
        if (is_small()) return Support(0);

        auto words = std::make_unique<uint64_t[]>(word_count_);
        const auto address = reinterpret_cast<std::uintptr_t>(words.get());
        allocations_.push_back(std::move(words));
        return Support(static_cast<uint64_t>(address));
    }

    Support clone(const Support& source)
    {
        Support result = make();
        copy(result, source);
        return result;
    }

    void copy(Support& destination, const Support& source) const noexcept
    {
        if (is_small()) destination.storage_ = source.storage_;
        else std::copy_n(words(source), word_count_, words(destination));
    }

    void swap(Support& left, Support& right) const noexcept { std::swap(left.storage_, right.storage_); }

    bool empty(const Support& value) const noexcept
    {
        if (is_small()) return value.storage_ == 0;
        for (size_t i = 0; i < word_count_; ++i)
            if (words(value)[i] != 0) return false;
        return true;
    }

    void clear(Support& value) const noexcept
    {
        if (is_small()) value.storage_ = 0;
        else std::fill_n(words(value), word_count_, uint64_t{0});
    }

    void set_all(Support& value) const noexcept
    {
        if (is_small()) {
            value.storage_ = set_all_n_bits(dimension_);
            return;
        }
        std::fill_n(words(value), word_count_, ~uint64_t{0});
        words(value)[word_count_ - 1] &= last_word_mask_;
    }

    bool contains(const Support& value, size_t position) const noexcept
    {
        assert(position < dimension_);
        if (is_small()) return (value.storage_ & (uint64_t{1} << position)) != 0;
        return (words(value)[position / 64] & (uint64_t{1} << (position % 64))) != 0;
    }

    void set(Support& value, size_t position) const noexcept
    {
        assert(position < dimension_);
        if (is_small()) value.storage_ |= uint64_t{1} << position;
        else words(value)[position / 64] |= uint64_t{1} << (position % 64);
    }

    void reset(Support& value, size_t position) const noexcept
    {
        assert(position < dimension_);
        if (is_small()) value.storage_ &= ~(uint64_t{1} << position);
        else words(value)[position / 64] &= ~(uint64_t{1} << (position % 64));
    }

    void add(Support& destination, const Support& source) const noexcept
    {
        if (is_small()) {
            destination.storage_ |= source.storage_;
            return;
        }
        for (size_t i = 0; i < word_count_; ++i) words(destination)[i] |= words(source)[i];
    }

    void subtract(Support& destination, const Support& source) const noexcept
    {
        if (is_small()) {
            destination.storage_ &= ~source.storage_;
            return;
        }
        for (size_t i = 0; i < word_count_; ++i) words(destination)[i] &= ~words(source)[i];
    }

    bool is_subset_of(const Support& subset, const Support& superset) const noexcept
    {
        if (is_small()) return (subset.storage_ & ~superset.storage_) == 0;
        for (size_t i = 0; i < word_count_; ++i)
            if ((words(subset)[i] & ~words(superset)[i]) != 0) return false;
        return true;
    }

    bool equal(const Support& left, const Support& right) const noexcept
    {
        if (is_small()) return left.storage_ == right.storage_;
        return std::equal(words(left), words(left) + word_count_, words(right));
    }

    bool less(const Support& left, const Support& right) const noexcept
    {
        if (is_small()) return left.storage_ < right.storage_;
        for (size_t i = word_count_; i > 0; --i) {
            if (words(left)[i - 1] != words(right)[i - 1]) return words(left)[i - 1] < words(right)[i - 1];
        }
        return false;
    }

    size_t count(const Support& value) const noexcept
    {
        if (is_small()) return popcount64(value.storage_);
        size_t result = 0;
        for (size_t i = 0; i < word_count_; ++i) result += popcount64(words(value)[i]);
        return result;
    }

    // Caller precondition: value is nonempty.
    size_t first(const Support& value) const noexcept
    {
        if (is_small()) return ctz64(value.storage_);
        for (size_t i = 0; i < word_count_; ++i)
            if (words(value)[i] != 0) return i * 64 + ctz64(words(value)[i]);
        assert(false && "SupportContext::first requires a nonempty support");
        return dimension_;
    }

    void extract_set_indices(const Support& value, std::vector<size_t>& indices) const
    {
        indices.clear();
        if (is_small()) {
            uint64_t word = value.storage_;
            while (word != 0) {
                indices.push_back(ctz64(word));
                word &= word - 1;
            }
            return;
        }

        for (size_t word_index = 0; word_index < word_count_; ++word_index) {
            uint64_t word = words(value)[word_index];
            while (word != 0) {
                indices.push_back(word_index * 64 + ctz64(word));
                word &= word - 1;
            }
        }
    }

    void extract_unset_indices(const Support& value, std::vector<size_t>& indices) const
    {
        indices.clear();
        if (is_small()) {
            uint64_t word = ~value.storage_ & set_all_n_bits(dimension_);
            while (word != 0) {
                indices.push_back(ctz64(word));
                word &= word - 1;
            }
            return;
        }

        for (size_t word_index = 0; word_index < word_count_; ++word_index) {
            uint64_t word = ~words(value)[word_index];
            if (word_index + 1 == word_count_) word &= last_word_mask_;
            while (word != 0) {
                indices.push_back(word_index * 64 + ctz64(word));
                word &= word - 1;
            }
        }
    }

    void rotate_one_right(Support& value) const noexcept
    {
        if (is_small()) {
            value.storage_ = rot_one_right(value.storage_, dimension_);
            return;
        }

        uint64_t* value_words = words(value);
        const bool wrap = (value_words[0] & uint64_t{1}) != 0;
        for (size_t i = 0; i + 1 < word_count_; ++i)
            value_words[i] = (value_words[i] >> 1) | (value_words[i + 1] << 63);
        value_words[word_count_ - 1] >>= 1;
        if (wrap) value_words[(dimension_ - 1) / 64] |= uint64_t{1} << ((dimension_ - 1) % 64);
    }

    void reflect(Support& value) const noexcept
    {
        if (is_small()) {
            value.storage_ = fracessa::support::reflect(value.storage_, dimension_);
            return;
        }

        uint64_t* value_words = words(value);
        size_t left = 0;
        size_t right = word_count_ - 1;
        while (left < right) {
            const uint64_t low = reverse_bits(value_words[left]);
            value_words[left++] = reverse_bits(value_words[right]);
            value_words[right--] = low;
        }
        if (left == right) value_words[left] = reverse_bits(value_words[left]);

        const size_t used_bits = dimension_ % 64;
        if (used_bits == 0) return;

        const size_t unused_bits = 64 - used_bits;
        for (size_t i = 0; i + 1 < word_count_; ++i)
            value_words[i] = (value_words[i] >> unused_bits) | (value_words[i + 1] << used_bits);
        value_words[word_count_ - 1] >>= unused_bits;
    }

    std::string to_bitstring(const Support& value) const
    {
        std::string result;
        result.reserve(dimension_);
        for (size_t position = dimension_; position > 0; --position)
            result += contains(value, position - 1) ? '1' : '0';
        return result;
    }

    std::string to_index_set(const Support& value) const
    {
        std::string result = "{";
        bool first_index = true;
        for (size_t position = 0; position < dimension_; ++position) {
            if (!contains(value, position)) continue;
            if (!first_index) result += ", ";
            result += std::to_string(position);
            first_index = false;
        }
        result += '}';
        return result;
    }

private:
    static uint64_t make_last_word_mask(size_t dimension) noexcept
    {
        const size_t used_bits = dimension % 64;
        return used_bits == 0 ? ~uint64_t{0} : (uint64_t{1} << used_bits) - 1;
    }

    static uint64_t* words(Support& value) noexcept
    {
        return reinterpret_cast<uint64_t*>(static_cast<std::uintptr_t>(value.storage_));
    }

    static const uint64_t* words(const Support& value) noexcept
    {
        return reinterpret_cast<const uint64_t*>(static_cast<std::uintptr_t>(value.storage_));
    }

    size_t dimension_;
    size_t word_count_;
    uint64_t last_word_mask_;
    std::vector<std::unique_ptr<uint64_t[]>> allocations_;
};

} // namespace fracessa::support
