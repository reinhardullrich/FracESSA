#pragma once

#include <array>
#include <cstddef>
#include <numeric>

#include <fracessa/bitset64.hpp>
#include <linalg/matrix_fraction.hpp>

/*
 * Exact extra symmetries of a symmetric circulant game.
 *
 * Bracelets already identify rotations and reflections. A unit a modulo n can
 * identify additional bracelets when multiplying every strategy index by a
 * leaves the complete matrix unchanged. The class detects those multipliers
 * once and filters only at the bracelet callback boundary.
 */
class CircularAffineSymmetry {
private:
    static constexpr size_t kMaxMultiplierClasses = (bs64::kMaxBitsetDimension + 1) / 2;

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

    inline bitset64 canonical_dihedral(bitset64 support) const noexcept {
        bitset64 smallest = support;
        bitset64 current = support;
        for (size_t shift = 1; shift < dimension_; ++shift) {
            current = bs64::rot_one_right(current, dimension_);
            if (current < smallest) smallest = current;
        }

        current = bs64::reflect(support, dimension_);
        for (size_t shift = 0; shift < dimension_; ++shift) {
            if (current < smallest) smallest = current;
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

    // V3 already emitted the smallest member of the identity multiplier's dihedral orbit. Reject as soon as another verified
    // multiplier reaches a smaller concrete support.
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
            const bitset64 image = canonical_dihedral(transform(support, multiplier_class));
            bool seen = false;
            for (size_t i = 0; i < image_count; ++i) {
                if (images[i] == image) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                images[image_count++] = image;
                callback(image);
            }
        }
    }
};
