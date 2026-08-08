#pragma once

#include <array>
#include <cstddef>
#include <numeric>
#include <vector>

#include <fracessa/bitset64.hpp>
#include <fracessa/bitset_multiword.hpp>
#include <linalg/matrix_fraction.hpp>

/*
 * Exact extra symmetries of a symmetric circulant game.
 *
 * Bracelets already identify rotations and reflections. A unit a modulo n is an integer coprime to n; it relabels strategy i as
 * a*i modulo n. When that relabeling leaves every matrix entry unchanged, it identifies additional equivalent bracelets. This class
 * detects those matrix-preserving multipliers once and filters only at the bracelet callback boundary.
 */
class CircularAffineSymmetry {
private:
    static constexpr size_t kMaxMultiplierClasses = (bs64::kMaxBitsetDimension + 1) / 2;

    struct CanonicalImage {
        bitset64 support;
        size_t right_shifts;
        bool reflected;
    };

    size_t dimension_;
    size_t multiplier_class_count_ = 0;
    std::array<std::array<bitset64, bs64::kMaxBitsetDimension>, kMaxMultiplierClasses> destination_bits_;

    inline void add_multiplier(size_t multiplier) noexcept {
        auto& destinations = destination_bits_[multiplier_class_count_++];
        for (size_t strategy = 0; strategy < dimension_; ++strategy)
            destinations[strategy] = 1ULL << ((multiplier * strategy) % dimension_);
    }

    inline bitset64 transform(bitset64 support, size_t multiplier_class) const noexcept {
        bitset64 transformed = 0;
        while (support) {
            transformed |= destination_bits_[multiplier_class][ctz64(support)];
            support &= support - 1;
        }
        return transformed;
    }

    inline CanonicalImage canonical_dihedral(bitset64 support) const noexcept {
        CanonicalImage smallest{support, 0, false};
        bitset64 current = support;
        for (size_t shift = 1; shift < dimension_; ++shift) {
            current = bs64::rot_one_right(current, dimension_);
            if (current < smallest.support) smallest = {current, shift, false};
        }

        current = bs64::reflect(support, dimension_);
        for (size_t shift = 0; shift < dimension_; ++shift) {
            if (current < smallest.support) smallest = {current, shift, true};
            current = bs64::rot_one_right(current, dimension_);
        }
        return smallest;
    }

public:
    explicit CircularAffineSymmetry(const linalg::matrix_frc& matrix) : dimension_(matrix.rows()) {
        add_multiplier(1);
        for (size_t multiplier = 2; multiplier <= dimension_ / 2; ++multiplier) {
            if (std::gcd(multiplier, dimension_) != 1) continue;

            bool preserves_matrix = true;
            for (size_t offset = 0; offset < dimension_; ++offset) {
                if (!(matrix(0, offset) == matrix(0, (multiplier * offset) % dimension_))) {
                    preserves_matrix = false;
                    break;
                }
            }
            if (preserves_matrix) add_multiplier(multiplier);
        }
    }

    inline size_t multiplier_class_count() const noexcept {
        return multiplier_class_count_;
    }

    inline size_t image_position(size_t position, size_t multiplier_class, bool reflected, size_t right_shifts) const noexcept {
        size_t image = ctz64(destination_bits_[multiplier_class][position]);
        if (reflected) image = dimension_ - 1 - image;
        return (image + dimension_ - right_shifts) % dimension_;
    }

    inline bitset64 image_mask(bitset64 support, size_t multiplier_class, bool reflected, size_t right_shifts) const noexcept {
        bitset64 image = transform(support, multiplier_class);
        if (reflected) image = bs64::reflect(image, dimension_);
        for (size_t shift = 0; shift < right_shifts; ++shift)
            image = bs64::rot_one_right(image, dimension_);
        return image;
    }

    // The bracelet generator already emitted the smallest rotation or reflection for the identity multiplier. Reject this support
    // if an additional matrix-preserving multiplier reaches a smaller concrete support.
    inline bool is_representative(bitset64 support) const noexcept {
        for (size_t multiplier_class = 1; multiplier_class < multiplier_class_count_; ++multiplier_class) {
            const bitset64 image = transform(support, multiplier_class);
            bitset64 transformed = image;
            for (size_t shift = 0; shift < dimension_; ++shift) {
                if (transformed < support) return false;
                transformed = bs64::rot_one_right(transformed, dimension_);
            }

            transformed = bs64::reflect(image, dimension_);
            for (size_t shift = 0; shift < dimension_; ++shift) {
                if (transformed < support) return false;
                transformed = bs64::rot_one_right(transformed, dimension_);
            }
        }
        return true;
    }

    // Return one representative for each distinct dihedral orbit inside the complete detected affine orbit.
    template<class Callback>
    inline void for_each_distinct_bracelet_image(bitset64 support, Callback&& callback) const {
        std::array<bitset64, kMaxMultiplierClasses> images{};
        size_t image_count = 0;

        for (size_t multiplier_class = 0; multiplier_class < multiplier_class_count_; ++multiplier_class) {
            const CanonicalImage image = canonical_dihedral(transform(support, multiplier_class));
            bool seen = false;
            for (size_t i = 0; i < image_count; ++i) {
                if (images[i] == image.support) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                images[image_count++] = image.support;
                callback(image.support, multiplier_class, image.reflected, image.right_shifts);
            }
        }
    }
};

/*
 * Multiword affine symmetry for the internal large circular path. Destination positions and all temporary masks are allocated once
 * in the constructor, so filtering a generated bracelet does not allocate.
 */
class CircularAffineSymmetryMultiword {
private:
    struct CanonicalImage {
        size_t right_shifts;
        bool reflected;
    };

    size_t dimension_;
    std::vector<std::vector<size_t>> destination_positions_;
    mutable std::vector<size_t> source_positions_;
    mutable bitset_multiword transformed_;
    mutable bitset_multiword canonical_;
    mutable std::vector<bitset_multiword> images_;

    inline void add_multiplier(size_t multiplier)
    {
        destination_positions_.emplace_back(dimension_);
        std::vector<size_t>& destinations = destination_positions_.back();
        for (size_t strategy = 0; strategy < dimension_; ++strategy)
            destinations[strategy] = (multiplier * strategy) % dimension_;
    }

    inline void transform_into(const bitset_multiword& support, size_t multiplier_class,
                               bitset_multiword& transformed) const
    {
        transformed.clear_all();
        support.extract_set_indices(source_positions_);
        const std::vector<size_t>& destinations = destination_positions_[multiplier_class];
        for (const size_t position : source_positions_) transformed.set_bit_at_pos(destinations[position]);
    }

    inline CanonicalImage canonical_dihedral(const bitset_multiword& support, size_t multiplier_class,
                                             bitset_multiword& smallest) const
    {
        CanonicalImage result{0, false};
        transform_into(support, multiplier_class, smallest);
        transformed_.copy_from(smallest);
        for (size_t shift = 1; shift < dimension_; ++shift) {
            transformed_.rotate_one_right();
            if (transformed_ < smallest) {
                smallest.copy_from(transformed_);
                result = {shift, false};
            }
        }

        // The loop stopped one rotation short of a complete cycle; this returns to the original transformed image.
        transformed_.rotate_one_right();
        transformed_.reflect();
        for (size_t shift = 0; shift < dimension_; ++shift) {
            if (transformed_ < smallest) {
                smallest.copy_from(transformed_);
                result = {shift, true};
            }
            transformed_.rotate_one_right();
        }
        return result;
    }

public:
    explicit CircularAffineSymmetryMultiword(const linalg::matrix_frc& matrix)
        : dimension_(matrix.rows()), transformed_(dimension_), canonical_(dimension_)
    {
        destination_positions_.reserve((dimension_ + 1) / 2);
        source_positions_.reserve(dimension_);

        add_multiplier(1);
        for (size_t multiplier = 2; multiplier <= dimension_ / 2; ++multiplier) {
            if (std::gcd(multiplier, dimension_) != 1) continue;

            bool preserves_matrix = true;
            for (size_t offset = 0; offset < dimension_; ++offset) {
                if (!(matrix(0, offset) == matrix(0, (multiplier * offset) % dimension_))) {
                    preserves_matrix = false;
                    break;
                }
            }
            if (preserves_matrix) add_multiplier(multiplier);
        }

        images_.reserve(destination_positions_.size());
        for (size_t i = 0; i < destination_positions_.size(); ++i) images_.emplace_back(dimension_);
    }

    inline size_t multiplier_class_count() const noexcept
    {
        return destination_positions_.size();
    }

    inline size_t image_position(size_t position, size_t multiplier_class, bool reflected, size_t right_shifts) const noexcept
    {
        size_t image = destination_positions_[multiplier_class][position];
        if (reflected) image = dimension_ - 1 - image;
        return (image + dimension_ - right_shifts) % dimension_;
    }

    inline void image_mask(const bitset_multiword& support, size_t multiplier_class, bool reflected,
                           size_t right_shifts, bitset_multiword& image) const
    {
        transform_into(support, multiplier_class, transformed_);
        if (reflected) transformed_.reflect();
        for (size_t shift = 0; shift < right_shifts; ++shift) transformed_.rotate_one_right();
        image.copy_from(transformed_);
    }

    // The bracelet generator already emitted the smallest rotation or reflection for the identity multiplier. Check only the
    // additional matrix-preserving multipliers for a smaller support.
    inline bool is_representative(const bitset_multiword& support) const
    {
        for (size_t multiplier_class = 1; multiplier_class < destination_positions_.size(); ++multiplier_class) {
            transform_into(support, multiplier_class, transformed_);
            for (size_t shift = 0; shift < dimension_; ++shift) {
                if (transformed_ < support) return false;
                transformed_.rotate_one_right();
            }

            // A complete rotation cycle restored the transformed image.
            transformed_.reflect();
            for (size_t shift = 0; shift < dimension_; ++shift) {
                if (transformed_ < support) return false;
                transformed_.rotate_one_right();
            }
        }
        return true;
    }

    // Callback references point into reusable class-owned storage and are valid only until this function is called again.
    template<class Callback>
    inline void for_each_distinct_bracelet_image(const bitset_multiword& support, Callback&& callback) const
    {
        size_t image_count = 0;
        for (size_t multiplier_class = 0; multiplier_class < destination_positions_.size(); ++multiplier_class) {
            const CanonicalImage image = canonical_dihedral(support, multiplier_class, canonical_);
            bool seen = false;
            for (size_t i = 0; i < image_count; ++i) {
                if (images_[i] == canonical_) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                images_[image_count].copy_from(canonical_);
                callback(images_[image_count++], multiplier_class, image.reflected, image.right_shifts);
            }
        }
    }
};
