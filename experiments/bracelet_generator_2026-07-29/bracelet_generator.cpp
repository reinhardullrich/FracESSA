// Standalone experiment: generate one circular support representative per
// rotation-and-reflection orbit (a fixed-weight binary bracelet).

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

class BraceletGenerator {
public:
    explicit BraceletGenerator(size_t dimension)
        : dimension_(dimension), word_(dimension + 1, 0)
    {}

    const std::vector<uint64_t>& generate(size_t support_size) {
        support_size_ = support_size;
        necklace_count_ = 0;
        bracelets_.clear();
        std::fill(word_.begin(), word_.end(), 0);
        generate_necklaces(1, 1, 0);
        std::sort(bracelets_.begin(), bracelets_.end());
        return bracelets_;
    }

    size_t necklace_count() const noexcept {
        return necklace_count_;
    }

    size_t orbit_size(uint64_t support) const noexcept {
        size_t rotations = 1;
        uint64_t rotated = rotate_right(support);
        while (rotated != support) {
            ++rotations;
            rotated = rotate_right(rotated);
        }

        const uint64_t reflected = reflect(support);
        rotated = support;
        for (size_t i = 0; i < rotations; ++i) {
            if (rotated == reflected)
                return rotations;
            rotated = rotate_right(rotated);
        }
        return 2 * rotations;
    }

private:
    size_t dimension_;
    size_t support_size_ = 0;
    size_t necklace_count_ = 0;
    std::vector<uint8_t> word_;
    std::vector<uint64_t> bracelets_;

    // Fredricksen-Kessler-Maiorana recursion emits one binary necklace per
    // rotational orbit. Tracking the weight avoids generating other sizes.
    void generate_necklaces(size_t position, size_t period, size_t ones) {
        if (position > dimension_) {
            if (dimension_ % period == 0 && ones == support_size_)
                emit_necklace();
            return;
        }
        if (ones > support_size_
            || ones + dimension_ - position + 1 < support_size_)
            return;

        word_[position] = word_[position - period];
        generate_necklaces(position + 1, period, ones + word_[position]);

        if (word_[position - period] == 0) {
            word_[position] = 1;
            generate_necklaces(position + 1, position, ones + 1);
        }
    }

    void emit_necklace() {
        ++necklace_count_;
        uint64_t support = 0;
        for (size_t i = 1; i <= dimension_; ++i) {
            if (word_[i] != 0)
                support |= 1ULL << (dimension_ - i);
        }

        const uint64_t rotation_rep = smallest_rotation(support);
        const uint64_t reflected_rep = smallest_rotation(reflect(support));
        if (rotation_rep <= reflected_rep)
            bracelets_.push_back(rotation_rep);
    }

    uint64_t rotate_right(uint64_t support) const noexcept {
        return (support >> 1) | ((support & 1ULL) << (dimension_ - 1));
    }

    uint64_t smallest_rotation(uint64_t support) const noexcept {
        uint64_t smallest = support;
        for (size_t i = 1; i < dimension_; ++i) {
            support = rotate_right(support);
            smallest = std::min(smallest, support);
        }
        return smallest;
    }

    uint64_t reflect(uint64_t support) const noexcept {
        uint64_t reflected = 0;
        for (size_t i = 0; i < dimension_; ++i)
            reflected |= ((support >> i) & 1ULL) << (dimension_ - 1 - i);
        return reflected;
    }
};

static uint64_t binomial_coefficient(size_t n, size_t k) noexcept {
    k = std::min(k, n - k);
    uint64_t result = 1;
    for (size_t i = 1; i <= k; ++i) {
        uint64_t numerator = static_cast<uint64_t>(n - k + i);
        uint64_t denominator = static_cast<uint64_t>(i);
        const uint64_t divisor = std::gcd(numerator, denominator);
        numerator /= divisor;
        denominator /= divisor;
        result /= denominator; // The remaining denominator divides C(n-k+i-1,i-1).
        result *= numerator;
    }
    return result;
}

static std::string bit_string(uint64_t support, size_t dimension) {
    std::string result;
    result.reserve(dimension);
    for (size_t i = dimension; i-- > 0;)
        result += ((support >> i) & 1ULL) != 0 ? '1' : '0';
    return result;
}

static void print_indices(uint64_t support) {
    std::cout << '{';
    bool first = true;
    while (support != 0) {
        const size_t index = static_cast<size_t>(__builtin_ctzll(support));
        if (!first)
            std::cout << ',';
        std::cout << index;
        first = false;
        support &= support - 1;
    }
    std::cout << '}';
}

static bool print_dimension(size_t dimension) {
    BraceletGenerator generator(dimension);
    uint64_t represented_total = 0;
    size_t bracelet_total = 0;

    std::cout << "dimension=" << dimension
              << " bit_order=" << dimension - 1 << "..0\n";
    for (size_t support_size = 1; support_size <= dimension; ++support_size) {
        const auto& bracelets = generator.generate(support_size);
        const uint64_t expected = binomial_coefficient(dimension, support_size);
        uint64_t represented = 0;
        for (const uint64_t support : bracelets)
            represented += generator.orbit_size(support);

        std::cout << "\nsupport_size=" << support_size
                  << " raw=" << expected
                  << " necklaces=" << generator.necklace_count()
                  << " bracelets=" << bracelets.size() << '\n';
        for (const uint64_t support : bracelets) {
            std::cout << "  mask=" << support
                      << " bits=" << bit_string(support, dimension)
                      << " indices=";
            print_indices(support);
            std::cout << " orbit=" << generator.orbit_size(support) << '\n';
        }

        if (represented != expected) {
            std::cerr << "orbit coverage failed for n=" << dimension
                      << ", support_size=" << support_size << '\n';
            return false;
        }
        represented_total += represented;
        bracelet_total += bracelets.size();
    }

    const uint64_t expected_total = (1ULL << dimension) - 1ULL;
    std::cout << "\nsummary bracelets=" << bracelet_total
              << " represented_nonempty_supports=" << represented_total
              << " expected=" << expected_total << '\n';
    return represented_total == expected_total;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " DIMENSION [DIMENSION ...]\n";
        return 2;
    }

    for (int i = 1; i < argc; ++i) {
        size_t dimension = 0;
        const std::string argument = argv[i];
        const auto parsed = std::from_chars(
            argument.data(), argument.data() + argument.size(), dimension);
        if (parsed.ec != std::errc{} || parsed.ptr != argument.data() + argument.size()
            || dimension == 0 || dimension >= 64) {
            std::cerr << "Dimension must be an integer from 1 through 63: "
                      << argument << '\n';
            return 2;
        }

        if (i > 1)
            std::cout << "\n============================================================\n\n";
        if (!print_dimension(dimension))
            return 1;
    }
    return 0;
}
