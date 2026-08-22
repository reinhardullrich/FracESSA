#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include <fracessa/support.hpp>
#include <fracessa/types.hpp>

namespace fracessa::support {

/** Exact extra symmetries of a symmetric circulant game. */
class CircularAffineSymmetry {
private:
    static constexpr size_t kMaxMultiplierClasses = (kMaxBitsetDimension + 1) / 2;

    struct SmallCanonicalImage {
        uint64_t support;
        size_t right_shifts;
        bool reflected;
    };

    struct LargeCanonicalImage {
        size_t right_shifts;
        bool reflected;
    };

    SupportContext& context_;
    size_t dimension_;
    size_t multiplier_class_count_ = 0;
    std::array<std::array<uint64_t, kMaxBitsetDimension>, kMaxMultiplierClasses> destination_bits_;
    std::vector<std::vector<size_t>> destination_positions_;
    mutable std::vector<size_t> source_positions_;
    mutable Support transformed_;
    mutable Support canonical_;
    mutable Support source_;
    mutable Support image_;
    mutable std::vector<Support> images_;

    inline void add_multiplier(size_t multiplier)
    {
        if (context_.is_small()) {
            auto& destinations = destination_bits_[multiplier_class_count_];
            for (size_t strategy = 0; strategy < dimension_; ++strategy)
                destinations[strategy] = uint64_t{1} << ((multiplier * strategy) % dimension_);
        } else {
            destination_positions_.emplace_back(dimension_);
            std::vector<size_t>& destinations = destination_positions_.back();
            for (size_t strategy = 0; strategy < dimension_; ++strategy)
                destinations[strategy] = (multiplier * strategy) % dimension_;
        }
        ++multiplier_class_count_;
    }

    inline uint64_t transform_small(uint64_t support, size_t multiplier_class) const noexcept
    {
        uint64_t transformed = 0;
        while (support != 0) {
            transformed |= destination_bits_[multiplier_class][ctz64(support)];
            support &= support - 1;
        }
        return transformed;
    }

    inline SmallCanonicalImage canonical_dihedral_small(uint64_t support) const noexcept
    {
        SmallCanonicalImage smallest{support, 0, false};
        uint64_t current = support;
        for (size_t shift = 1; shift < dimension_; ++shift) {
            current = rot_one_right(current, dimension_);
            if (current < smallest.support) smallest = {current, shift, false};
        }

        current = reflect(support, dimension_);
        for (size_t shift = 0; shift < dimension_; ++shift) {
            if (current < smallest.support) smallest = {current, shift, true};
            current = rot_one_right(current, dimension_);
        }
        return smallest;
    }

    inline void transform_large(const Support& support, size_t multiplier_class, Support& transformed) const
    {
        context_.clear(transformed);
        context_.extract_set_indices(support, source_positions_);
        const std::vector<size_t>& destinations = destination_positions_[multiplier_class];
        for (const size_t position : source_positions_) context_.set(transformed, destinations[position]);
    }

    inline LargeCanonicalImage canonical_dihedral_large(const Support& support, size_t multiplier_class,
                                                         Support& smallest) const
    {
        LargeCanonicalImage result{0, false};
        transform_large(support, multiplier_class, smallest);
        context_.copy(transformed_, smallest);
        for (size_t shift = 1; shift < dimension_; ++shift) {
            context_.rotate_one_right(transformed_);
            if (context_.less(transformed_, smallest)) {
                context_.copy(smallest, transformed_);
                result = {shift, false};
            }
        }

        context_.rotate_one_right(transformed_);
        context_.reflect(transformed_);
        for (size_t shift = 0; shift < dimension_; ++shift) {
            if (context_.less(transformed_, smallest)) {
                context_.copy(smallest, transformed_);
                result = {shift, true};
            }
            context_.rotate_one_right(transformed_);
        }
        return result;
    }

public:
    CircularAffineSymmetry(const numeric::matrix_int& matrix, SupportContext& context)
        : context_(context)
        , dimension_(context.dimension())
        , transformed_(context.make())
        , canonical_(context.make())
        , source_(context.make())
        , image_(context.make())
    {
        if (!context_.is_small()) {
            destination_positions_.reserve((dimension_ + 1) / 2);
            source_positions_.reserve(dimension_);
        }

        add_multiplier(1);
        for (size_t multiplier = 2; multiplier <= dimension_ / 2; ++multiplier) {
            if (std::gcd(multiplier, dimension_) != 1) continue;

            bool preserves_matrix = true;
            for (size_t offset = 0; offset < dimension_; ++offset) {
                if (matrix(0, offset).compare(matrix(0, (multiplier * offset) % dimension_)) != 0) {
                    preserves_matrix = false;
                    break;
                }
            }
            if (preserves_matrix) add_multiplier(multiplier);
        }

        if (!context_.is_small()) {
            images_.reserve(multiplier_class_count_);
            for (size_t i = 0; i < multiplier_class_count_; ++i) images_.push_back(context_.make());
        }
    }

    inline size_t multiplier_class_count() const noexcept { return multiplier_class_count_; }

    inline size_t image_position(size_t position, size_t multiplier_class, bool reflected, size_t right_shifts) const noexcept
    {
        size_t image = context_.is_small() ? ctz64(destination_bits_[multiplier_class][position])
                                           : destination_positions_[multiplier_class][position];
        if (reflected) image = dimension_ - 1 - image;
        return (image + dimension_ - right_shifts) % dimension_;
    }

    inline void image_mask(const Support& support, size_t multiplier_class, bool reflected,
                           size_t right_shifts, Support& image) const
    {
        if (context_.is_small()) {
            uint64_t bits = transform_small(context_.small_bits(support), multiplier_class);
            if (reflected) bits = reflect(bits, dimension_);
            for (size_t shift = 0; shift < right_shifts; ++shift) bits = rot_one_right(bits, dimension_);
            context_.set_small_bits(image, bits);
            return;
        }

        transform_large(support, multiplier_class, transformed_);
        if (reflected) context_.reflect(transformed_);
        for (size_t shift = 0; shift < right_shifts; ++shift) context_.rotate_one_right(transformed_);
        context_.copy(image, transformed_);
    }

    inline bool is_representative(const Support& support) const
    {
        if (context_.is_small()) {
            const uint64_t support_bits = context_.small_bits(support);
            for (size_t multiplier_class = 1; multiplier_class < multiplier_class_count_; ++multiplier_class) {
                const uint64_t image = transform_small(support_bits, multiplier_class);
                uint64_t transformed = image;
                for (size_t shift = 0; shift < dimension_; ++shift) {
                    if (transformed < support_bits) return false;
                    transformed = rot_one_right(transformed, dimension_);
                }

                transformed = reflect(image, dimension_);
                for (size_t shift = 0; shift < dimension_; ++shift) {
                    if (transformed < support_bits) return false;
                    transformed = rot_one_right(transformed, dimension_);
                }
            }
            return true;
        }

        for (size_t multiplier_class = 1; multiplier_class < multiplier_class_count_; ++multiplier_class) {
            transform_large(support, multiplier_class, transformed_);
            for (size_t shift = 0; shift < dimension_; ++shift) {
                if (context_.less(transformed_, support)) return false;
                context_.rotate_one_right(transformed_);
            }

            context_.reflect(transformed_);
            for (size_t shift = 0; shift < dimension_; ++shift) {
                if (context_.less(transformed_, support)) return false;
                context_.rotate_one_right(transformed_);
            }
        }
        return true;
    }

    template<class Callback>
    inline void for_each_distinct_bracelet_image(const Support& support, Callback&& callback) const
    {
        if (context_.is_small()) {
            const uint64_t source = context_.small_bits(support);
            std::array<uint64_t, kMaxMultiplierClasses> images{};
            size_t image_count = 0;
            for (size_t multiplier_class = 0; multiplier_class < multiplier_class_count_; ++multiplier_class) {
                const SmallCanonicalImage image = canonical_dihedral_small(transform_small(source, multiplier_class));
                bool seen = false;
                for (size_t i = 0; i < image_count; ++i) {
                    if (images[i] == image.support) {
                        seen = true;
                        break;
                    }
                }
                if (!seen) {
                    images[image_count++] = image.support;
                    context_.set_small_bits(image_, image.support);
                    callback(image_, multiplier_class, image.reflected, image.right_shifts);
                }
            }
            return;
        }

        context_.copy(source_, support);
        size_t image_count = 0;
        for (size_t multiplier_class = 0; multiplier_class < multiplier_class_count_; ++multiplier_class) {
            const LargeCanonicalImage image = canonical_dihedral_large(source_, multiplier_class, canonical_);
            bool seen = false;
            for (size_t i = 0; i < image_count; ++i) {
                if (context_.equal(images_[i], canonical_)) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                context_.copy(images_[image_count], canonical_);
                callback(images_[image_count++], multiplier_class, image.reflected, image.right_shifts);
            }
        }
    }
};

} // namespace fracessa::support
