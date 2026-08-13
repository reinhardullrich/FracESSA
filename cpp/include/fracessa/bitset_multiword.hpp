#pragma once

#include <fracessa/bitset.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace fracessa::support {

// Runtime-sized support storage. Production uses this type only above dimension 64;
// lower dimensions remain raw uint64_t so their hot path is unchanged.
class bitset_multiword {
public:
    explicit bitset_multiword(size_t dimension)
        : dimension_(dimension), words_(word_count_for(dimension), 0), last_word_mask_(make_last_word_mask(dimension))
    {
        if (dimension == 0) throw std::invalid_argument("bitset_multiword dimension must be positive");
    }

    size_t dimension() const noexcept { return dimension_; }
    size_t word_count() const noexcept { return words_.size(); }

    bool empty() const noexcept
    {
        for (const uint64_t word : words_)
            if (word != 0) return false;
        return true;
    }

    void clear_all() noexcept { std::fill(words_.begin(), words_.end(), uint64_t{0}); }

    void set_all() noexcept
    {
        std::fill(words_.begin(), words_.end(), ~uint64_t{0});
        words_.back() &= last_word_mask_;
    }

    void copy_from(const bitset_multiword& other) noexcept
    {
        assert_same_dimension(other);
        std::copy(other.words_.begin(), other.words_.end(), words_.begin());
    }

    bool is_set_at_pos(size_t position) const noexcept
    {
        assert(position < dimension_);
        return ((words_[position / 64] >> (position % 64)) & uint64_t{1}) != 0;
    }

    void set_bit_at_pos(size_t position) noexcept
    {
        assert(position < dimension_);
        words_[position / 64] |= uint64_t{1} << (position % 64);
    }

    void clear_bit_at_pos(size_t position) noexcept
    {
        assert(position < dimension_);
        words_[position / 64] &= ~(uint64_t{1} << (position % 64));
    }

    void union_with(const bitset_multiword& other) noexcept
    {
        assert_same_dimension(other);
        for (size_t i = 0; i < words_.size(); ++i) words_[i] |= other.words_[i];
    }

    void subtract(const bitset_multiword& other) noexcept
    {
        assert_same_dimension(other);
        for (size_t i = 0; i < words_.size(); ++i) words_[i] &= ~other.words_[i];
    }

    bool is_subset_of(const bitset_multiword& other) const noexcept
    {
        assert_same_dimension(other);
        for (size_t i = 0; i < words_.size(); ++i)
            if ((words_[i] & ~other.words_[i]) != 0) return false;
        return true;
    }

    size_t count_set_bits() const noexcept
    {
        size_t count = 0;
        for (const uint64_t word : words_) count += popcount64(word);
        return count;
    }

    // Caller precondition: the set is nonempty.
    size_t find_pos_first_set_bit() const noexcept
    {
        for (size_t i = 0; i < words_.size(); ++i)
            if (words_[i] != 0) return i * 64 + ctz64(words_[i]);
        assert(false && "find_pos_first_set_bit requires a nonempty set");
        return dimension_;
    }

    // The caller can reserve dimension() entries once and reuse the vector without allocations.
    void extract_set_indices(std::vector<size_t>& indices) const
    {
        indices.clear();
        for (size_t word_index = 0; word_index < words_.size(); ++word_index) {
            uint64_t word = words_[word_index];
            while (word != 0) {
                indices.push_back(word_index * 64 + ctz64(word));
                word &= word - 1;
            }
        }
    }

    // As above, but extract positions not present in the set without constructing its complement.
    void extract_unset_indices(std::vector<size_t>& indices) const
    {
        indices.clear();
        for (size_t word_index = 0; word_index < words_.size(); ++word_index) {
            uint64_t word = ~words_[word_index];
            if (word_index + 1 == words_.size()) word &= last_word_mask_;
            while (word != 0) {
                indices.push_back(word_index * 64 + ctz64(word));
                word &= word - 1;
            }
        }
    }

    // Move strategy i to i-1 modulo dimension, matching rot_one_right().
    void rotate_one_right() noexcept
    {
        const bool wrap = (words_[0] & uint64_t{1}) != 0;
        for (size_t i = 0; i + 1 < words_.size(); ++i)
            words_[i] = (words_[i] >> 1) | (words_[i + 1] << 63);
        words_.back() >>= 1;
        if (wrap) set_bit_at_pos(dimension_ - 1);
    }

    // Mirror strategy i to dimension-1-i one word at a time without allocating scratch storage.
    void reflect() noexcept
    {
        size_t left = 0;
        size_t right = words_.size() - 1;
        while (left < right) {
            const uint64_t low = reverse_bits(words_[left]);
            words_[left++] = reverse_bits(words_[right]);
            words_[right--] = low;
        }
        if (left == right) words_[left] = reverse_bits(words_[left]);

        const size_t used_bits = dimension_ % 64;
        if (used_bits == 0) return;

        const size_t unused_bits = 64 - used_bits;
        for (size_t i = 0; i + 1 < words_.size(); ++i)
            words_[i] = (words_[i] >> unused_bits) | (words_[i + 1] << used_bits);
        words_.back() >>= unused_bits;
    }

    std::string to_bitstring() const
    {
        std::string result;
        result.reserve(dimension_);
        for (size_t position = dimension_; position > 0; --position)
            result += is_set_at_pos(position - 1) ? '1' : '0';
        return result;
    }

    friend bool operator==(const bitset_multiword& left, const bitset_multiword& right) noexcept
    {
        return left.dimension_ == right.dimension_ && left.words_ == right.words_;
    }

    friend bool operator!=(const bitset_multiword& left, const bitset_multiword& right) noexcept { return !(left == right); }

    friend bool operator<(const bitset_multiword& left, const bitset_multiword& right) noexcept
    {
        left.assert_same_dimension(right);
        for (size_t i = left.words_.size(); i > 0; --i) {
            if (left.words_[i - 1] != right.words_[i - 1]) return left.words_[i - 1] < right.words_[i - 1];
        }
        return false;
    }

private:
    static size_t word_count_for(size_t dimension) noexcept
    {
        return dimension / 64 + static_cast<size_t>(dimension % 64 != 0);
    }

    static uint64_t make_last_word_mask(size_t dimension) noexcept
    {
        const size_t used_bits = dimension % 64;
        return used_bits == 0 ? ~uint64_t{0} : (uint64_t{1} << used_bits) - 1;
    }

    void assert_same_dimension(const bitset_multiword& other) const noexcept { assert(dimension_ == other.dimension_); }

    size_t dimension_;
    std::vector<uint64_t> words_;
    uint64_t last_word_mask_;
};

} // namespace fracessa::support
