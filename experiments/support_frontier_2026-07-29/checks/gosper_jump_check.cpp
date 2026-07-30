#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using mask = std::uint64_t;

size_t popcount(mask value) { return __builtin_popcountll(value); }
size_t ctz(mask value) { return __builtin_ctzll(value); }

mask next_same_size(mask value)
{
    const mask lowest = value & (0ULL - value);
    const mask ripple = value + lowest;
    return ripple | (((ripple ^ value) >> 2) >> ctz(lowest));
}

mask after_forbidden_block(mask value, size_t lowest, size_t dimension)
{
    const mask limit = 1ULL << dimension;
    const mask full = limit - 1;
    const mask lower = (1ULL << (lowest + 1)) - 1;
    const mask zeros_above = (~value & full) & ~lower;
    if (zeros_above == 0)
        return limit;

    const size_t pivot = ctz(zeros_above);
    const mask below_pivot = (1ULL << pivot) - 1;
    const size_t selected_below = popcount(value & below_pivot);
    return (value & ~((1ULL << (pivot + 1)) - 1))
        | (1ULL << pivot)
        | ((1ULL << (selected_below - 1)) - 1);
}

bool forbidden(mask value, const std::vector<mask>& forbidden_sets)
{
    for (const mask subset : forbidden_sets)
        if ((subset & ~value) == 0)
            return true;
    return false;
}

std::vector<mask> generate(size_t dimension, size_t cardinality,
                           const std::vector<mask>& forbidden_sets, bool jump)
{
    const mask limit = 1ULL << dimension;
    mask value = (1ULL << cardinality) - 1;
    std::vector<mask> result;

    while (value < limit) {
        if (!forbidden(value, forbidden_sets)) {
            result.push_back(value);
            value = next_same_size(value);
            continue;
        }

        if (!jump) {
            value = next_same_size(value);
            continue;
        }

        size_t best_lowest = 0;
        for (const mask subset : forbidden_sets)
            if ((subset & ~value) == 0 && ctz(subset) > best_lowest)
                best_lowest = ctz(subset);
        value = after_forbidden_block(value, best_lowest, dimension);
    }
    return result;
}

int main()
{
    // Every single forbidden subset through n=12 checks order and completeness.
    for (size_t n = 2; n <= 12; ++n) {
        const mask limit = 1ULL << n;
        for (mask subset = 1; subset < limit; ++subset) {
            for (size_t k = popcount(subset) + 1; k <= n; ++k) {
                const std::vector<mask> forbidden_sets{subset};
                assert(generate(n, k, forbidden_sets, false)
                    == generate(n, k, forbidden_sets, true));
            }
        }
    }

    // Overlapping forbidden families exercise witness selection and repeated jumps.
    std::mt19937_64 random(20260729);
    for (size_t trial = 0; trial < 2000; ++trial) {
        const size_t n = 6 + random() % 15;
        const size_t k = 2 + random() % (n - 1);
        std::vector<mask> forbidden_sets;
        for (size_t i = 0; i < 1 + random() % 20; ++i) {
            mask subset = random() & ((1ULL << n) - 1);
            if (subset != 0 && popcount(subset) < k)
                forbidden_sets.push_back(subset);
        }
        assert(generate(n, k, forbidden_sets, false)
            == generate(n, k, forbidden_sets, true));
    }

    std::cout << "gosper jump sequence matches filtered Gosper sequence\n";
}
