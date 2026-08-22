#include <fracessa/support_generator_sat.hpp>

#include <array>
#include <cassert>
#include <cstddef>

using fracessa::support::SatSupportGenerator;
using fracessa::support::Support;
using fracessa::support::SupportContext;

int main()
{
    {
        SupportContext context(4);
        SatSupportGenerator generator(context);
        std::array<size_t, 5> counts{};
        Support support = context.make();
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
        SupportContext context(4);
        SatSupportGenerator generator(context);
        Support required = context.make();
        Support invaders = context.make();
        context.set_small_bits(required, 0b0011);
        generator.add_upper_cone(required);
        context.set_small_bits(required, 0b0100);
        context.set_small_bits(invaders, 0b1000);
        generator.add_invader_interval(required, invaders);

        std::array<bool, 16> seen{};
        Support support = context.make();
        for (size_t cardinality = 1; cardinality <= 4; ++cardinality) {
            generator.start_cardinality(cardinality);
            while (generator.take_first(support)) {
                seen[context.small_bits(support)] = true;
                generator.add_exact(support, cardinality);
            }
        }
        for (size_t candidate = 1; candidate < seen.size(); ++candidate) {
            const bool in_upper_cone = (candidate & 0b0011) == 0b0011;
            const bool in_invader_interval = (candidate & 0b0100) != 0 && (candidate & 0b1000) == 0;
            assert(seen[candidate] == !(in_upper_cone || in_invader_interval));
        }
    }

    {
        SupportContext context(65);
        SatSupportGenerator generator(context);
        Support support = context.make();
        generator.start_cardinality(1);
        size_t count = 0;
        while (generator.take_first(support)) {
            assert(context.count(support) == 1);
            generator.add_exact(support, 1);
            ++count;
        }
        assert(count == 65);
    }
}
