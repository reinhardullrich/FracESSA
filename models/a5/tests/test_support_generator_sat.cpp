#include <fracessa/support_generator_sat.hpp>

#include <array>
#include <cassert>
#include <cstddef>

using fracessa::support::SatSupportGenerator;
using fracessa::support::bitset;
using fracessa::support::bitset_multiword;

int main()
{
    {
        SatSupportGenerator generator(4);
        std::array<size_t, 5> counts{};
        bitset support = 0;
        for (size_t cardinality = 1; cardinality <= 4; ++cardinality) {
            generator.start_cardinality(cardinality);
            while (generator.take_first(support)) {
                ++counts[cardinality];
                generator.add_exact(support, cardinality);
            }
        }
        assert((counts == std::array<size_t, 5>{0, 4, 6, 4, 1}));
    }

    {
        SatSupportGenerator generator(4);
        generator.add_upper_cone(bitset{0b0011});
        generator.add_invader_interval(bitset{0b0100}, bitset{0b1000});

        std::array<bool, 16> seen{};
        bitset support = 0;
        for (size_t cardinality = 1; cardinality <= 4; ++cardinality) {
            generator.start_cardinality(cardinality);
            while (generator.take_first(support)) {
                seen[support] = true;
                generator.add_exact(support, cardinality);
            }
        }
        for (bitset candidate = 1; candidate < seen.size(); ++candidate) {
            const bool in_upper_cone = (candidate & 0b0011) == 0b0011;
            const bool in_invader_interval = (candidate & 0b0100) != 0 && (candidate & 0b1000) == 0;
            assert(seen[candidate] == !(in_upper_cone || in_invader_interval));
        }
    }

    {
        SatSupportGenerator generator(65);
        bitset_multiword support(65);
        generator.start_cardinality(1);
        size_t count = 0;
        while (generator.take_first(support)) {
            assert(support.count_set_bits() == 1);
            generator.add_exact(support, 1);
            ++count;
        }
        assert(count == 65);
    }
}
