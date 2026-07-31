#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;

struct observation {
    uint64_t count = 0;
    uint64_t sum = 0;
};

struct measurement {
    observation value;
    uint64_t nanoseconds = 0;
};

using subset_generator = observation (*)(size_t);
using eager_buckets = std::vector<std::vector<uint64_t>>;

struct eager_profile {
    observation value;
    uint64_t creation_nanoseconds = 0;
    uint64_t iteration_nanoseconds = 0;
    uint64_t cleanup_nanoseconds = 0;
};

uint64_t binomial_coefficient(uint64_t n, uint64_t k)
{
    if (k > n) return 0;
    if (k > n - k) k = n - k;

    uint64_t result = 1;
    for (uint64_t i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

inline void observe(observation& result, uint64_t subset) noexcept
{
    ++result.count;
    result.sum += subset;
}

eager_buckets create_eager_buckets(size_t dimension)
{
    eager_buckets subsets(dimension);
    for (size_t cardinality = 1; cardinality <= dimension; ++cardinality) {
        subsets[cardinality - 1].reserve(
            binomial_coefficient(dimension, cardinality));
    }

    const uint64_t limit = uint64_t{1} << dimension;
    for (uint64_t subset = 1; subset < limit; ++subset) {
        const size_t cardinality =
            static_cast<size_t>(__builtin_popcountll(subset));
        subsets[cardinality - 1].push_back(subset);
    }

    return subsets;
}

observation iterate_eager_buckets(const eager_buckets& subsets)
{
    observation result;
    for (const auto& cardinality : subsets) {
        for (uint64_t subset : cardinality) observe(result, subset);
    }
    return result;
}

// Current copositivity behavior: build all cardinality buckets, then iterate.
observation eager(size_t dimension)
{
    const eager_buckets subsets = create_eager_buckets(dimension);
    return iterate_eager_buckets(subsets);
}

eager_profile profile_eager(size_t dimension)
{
    observation value;
    const auto start = clock_type::now();
    auto created = start;
    auto iterated = start;
    {
        const eager_buckets subsets = create_eager_buckets(dimension);
        created = clock_type::now();
        value = iterate_eager_buckets(subsets);
        iterated = clock_type::now();
    }
    const auto cleaned = clock_type::now();

    return {
        value,
        static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(created - start)
                .count()),
        static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(iterated - created)
                .count()),
        static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(cleaned - iterated)
                .count()),
    };
}

double median(std::vector<uint64_t>& values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return (static_cast<double>(values[middle - 1]) +
            static_cast<double>(values[middle])) / 2.0;
}

// Gosper's hack: stream the next integer with the same number of set bits.
observation gosper(size_t dimension)
{
    observation result;
    const uint64_t limit = uint64_t{1} << dimension;
    for (size_t cardinality = 1; cardinality <= dimension; ++cardinality) {
        uint64_t subset = (uint64_t{1} << cardinality) - 1;
        while (subset < limit) {
            observe(result, subset);
            const uint64_t lowest = subset & (uint64_t{0} - subset);
            const uint64_t ripple = subset + lowest;
            subset = ripple | (((ripple ^ subset) >> 2) / lowest);
        }
    }
    return result;
}

// Production-style bitwise DFS: exclude the current high bit, then include it.
void dfs_from(
    size_t bits_remaining,
    size_t needed,
    uint64_t partial,
    observation& result)
{
    if (needed == 0) {
        observe(result, partial);
        return;
    }
    if (needed > bits_remaining) return;

    const size_t bit = bits_remaining - 1;
    if (needed < bits_remaining) {
        dfs_from(bit, needed, partial, result);
    }
    dfs_from(bit, needed - 1, partial | (uint64_t{1} << bit), result);
}

observation dfs(size_t dimension)
{
    observation result;
    for (size_t cardinality = 1; cardinality <= dimension; ++cardinality) {
        dfs_from(dimension, cardinality, 0, result);
    }
    return result;
}

measurement measure(subset_generator generator, size_t dimension)
{
    const auto start = clock_type::now();
    observation value = generator(dimension);
    const auto end = clock_type::now();
    return {
        value,
        static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count())
    };
}

measurement measure_batch(
    subset_generator generator,
    size_t dimension,
    size_t repetitions)
{
    observation value;
    volatile size_t current_dimension = dimension;
    const auto start = clock_type::now();
    for (size_t i = 0; i < repetitions; ++i) {
        const observation current = generator(current_dimension);
        value.count += current.count;
        value.sum += current.sum;
    }
    const auto end = clock_type::now();
    return {
        value,
        static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count())
    };
}

void verify(size_t dimension, const std::array<observation, 3>& values)
{
    const uint64_t expected_count = (uint64_t{1} << dimension) - 1;
    const uint64_t expected_sum =
        expected_count * (uint64_t{1} << (dimension - 1));
    for (const observation& value : values) {
        if (value.count != expected_count || value.sum != expected_sum) {
            throw std::runtime_error("generator mismatch at dimension " +
                                     std::to_string(dimension));
        }
    }
}

void benchmark_median_range(size_t first_dimension, size_t last_dimension)
{
    constexpr size_t sample_count = 100;
    constexpr uint64_t target_subset_visits = 310'000;
    constexpr std::array<subset_generator, 3> generators{eager, gosper, dfs};

    std::cout
        << "dimension\tsamples\trepetitions_per_sample\tsubsets"
           "\teager_median_ns\tgosper_median_ns\tdfs_median_ns"
           "\tgosper_vs_eager_percent\tdfs_vs_eager_percent\tchecksum\n";

    for (size_t dimension = first_dimension;
         dimension <= last_dimension;
         ++dimension) {
        const uint64_t subsets = (uint64_t{1} << dimension) - 1;
        const size_t repetitions = static_cast<size_t>(
            (target_subset_visits + subsets - 1) / subsets);
        const uint64_t expected_count = subsets * repetitions;
        const uint64_t expected_sum =
            subsets * (uint64_t{1} << (dimension - 1)) * repetitions;
        const auto verify_batch = [expected_count, expected_sum](
                                      const measurement& measured) {
            if (measured.value.count != expected_count ||
                measured.value.sum != expected_sum) {
                throw std::runtime_error("generator batch mismatch");
            }
        };

        for (subset_generator generator : generators) {
            verify_batch(measure_batch(generator, dimension, repetitions));
        }

        std::array<std::vector<uint64_t>, 3> samples;
        for (auto& algorithm_samples : samples) {
            algorithm_samples.reserve(sample_count);
        }

        for (size_t sample = 0; sample < sample_count; ++sample) {
            std::array<size_t, 3> order{0, 1, 2};
            std::rotate(
                order.begin(),
                order.begin() + (sample + dimension) % order.size(),
                order.end());
            for (size_t algorithm : order) {
                const measurement measured =
                    measure_batch(generators[algorithm], dimension, repetitions);
                verify_batch(measured);
                samples[algorithm].push_back(measured.nanoseconds);
            }
        }

        std::array<double, 3> median_ns{};
        for (size_t algorithm = 0; algorithm < samples.size(); ++algorithm) {
            median_ns[algorithm] = median(samples[algorithm]) / repetitions;
        }

        std::cout
            << dimension << '\t'
            << sample_count << '\t'
            << repetitions << '\t'
            << subsets << '\t'
            << std::fixed << std::setprecision(6)
            << median_ns[0] << '\t'
            << median_ns[1] << '\t'
            << median_ns[2] << '\t'
            << 100.0 * (median_ns[1] / median_ns[0] - 1.0) << '\t'
            << 100.0 * (median_ns[2] / median_ns[0] - 1.0) << '\t'
            << expected_sum / repetitions << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::string_view(argv[1]) == "median-5") {
        benchmark_median_range(5, 5);
        return 0;
    }
    if (argc == 2 && std::string_view(argv[1]) == "median-3-25") {
        benchmark_median_range(3, 25);
        return 0;
    }
    if (argc == 2 && std::string_view(argv[1]) == "profile-eager-30") {
        const observation warm = eager(30);
        verify(30, {warm, warm, warm});

        const eager_profile profile = profile_eager(30);
        verify(30, {profile.value, profile.value, profile.value});
        const uint64_t total = profile.creation_nanoseconds +
                               profile.iteration_nanoseconds +
                               profile.cleanup_nanoseconds;
        const auto percentage = [total](uint64_t part) {
            return 100.0 * static_cast<double>(part) / static_cast<double>(total);
        };

        std::cout
            << "dimension\tsubsets\tcreation_ns\titeration_ns\tcleanup_ns\ttotal_ns"
               "\tcreation_percent\titeration_percent\tcleanup_percent\tchecksum\n"
            << 30 << '\t'
            << profile.value.count << '\t'
            << profile.creation_nanoseconds << '\t'
            << profile.iteration_nanoseconds << '\t'
            << profile.cleanup_nanoseconds << '\t'
            << total << '\t'
            << std::fixed << std::setprecision(6)
            << percentage(profile.creation_nanoseconds) << '\t'
            << percentage(profile.iteration_nanoseconds) << '\t'
            << percentage(profile.cleanup_nanoseconds) << '\t'
            << profile.value.sum << '\n';
        return 0;
    }
    if (argc != 1) {
        std::cerr
            << "Usage: subset_generation_benchmark "
               "[profile-eager-30|median-5|median-3-25]\n";
        return 2;
    }

    constexpr std::array<subset_generator, 3> generators{eager, gosper, dfs};

    std::cout
        << "dimension\tsubsets\teager_ns\tgosper_ns\tdfs_ns"
           "\tgosper_vs_eager_percent\tdfs_vs_eager_percent\tchecksum\n";

    for (size_t dimension = 2; dimension <= 30; ++dimension) {
        std::array<observation, 3> warm{};
        for (size_t i = 0; i < generators.size(); ++i) {
            warm[i] = generators[i](dimension);
        }
        verify(dimension, warm);

        std::array<size_t, 3> order{0, 1, 2};
        std::rotate(order.begin(), order.begin() + dimension % order.size(),
                    order.end());

        std::array<measurement, 3> measured{};
        for (size_t algorithm : order) {
            measured[algorithm] = measure(generators[algorithm], dimension);
        }
        verify(dimension, {
            measured[0].value,
            measured[1].value,
            measured[2].value,
        });

        const double eager_ns = static_cast<double>(measured[0].nanoseconds);
        const double gosper_change =
            100.0 * (static_cast<double>(measured[1].nanoseconds) / eager_ns - 1.0);
        const double dfs_change =
            100.0 * (static_cast<double>(measured[2].nanoseconds) / eager_ns - 1.0);

        std::cout << dimension << '\t'
                  << measured[0].value.count << '\t'
                  << measured[0].nanoseconds << '\t'
                  << measured[1].nanoseconds << '\t'
                  << measured[2].nanoseconds << '\t'
                  << std::fixed << std::setprecision(6)
                  << gosper_change << '\t'
                  << dfs_change << '\t'
                  << measured[0].value.sum << '\n';
    }

    std::cerr << "All three generators matched for dimensions 2-30.\n";
}
