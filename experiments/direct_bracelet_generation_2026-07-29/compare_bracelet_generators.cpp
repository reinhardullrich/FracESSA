#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

constexpr size_t kMaxDimension = 63;
using Clock = std::chrono::steady_clock;

uint64_t rotate_right(uint64_t support, size_t dimension) noexcept
{
    return (support >> 1) | ((support & 1ULL) << (dimension - 1));
}

uint64_t reflect(uint64_t support, size_t dimension) noexcept
{
    uint64_t reflected = 0;
    for (size_t i = 0; i < dimension; ++i) {
        reflected |= ((support >> i) & 1ULL) << (dimension - 1 - i);
    }
    return reflected;
}

uint64_t canonical_bracelet(uint64_t support, size_t dimension) noexcept
{
    uint64_t smallest = support;
    uint64_t mirrored = reflect(support, dimension);
    for (size_t i = 0; i < dimension; ++i) {
        smallest = std::min({ smallest, support, mirrored });
        support = rotate_right(support, dimension);
        mirrored = rotate_right(mirrored, dimension);
    }
    return smallest;
}

// Current experiment: generate fixed-weight necklaces, then discard the
// lexicographically larger of the two necklace classes related by reflection.
class FkmBraceletGenerator {
public:
    explicit FkmBraceletGenerator(size_t dimension)
        : dimension_(dimension)
        , word_(dimension + 1, 0)
    {}

    const std::vector<uint64_t>& generate(size_t support_size)
    {
        support_size_ = support_size;
        bracelets_.clear();
        std::fill(word_.begin(), word_.end(), 0);
        generate_necklaces(1, 1, 0);
        return bracelets_;
    }

private:
    size_t dimension_;
    size_t support_size_ = 0;
    std::vector<uint8_t> word_;
    std::vector<uint64_t> bracelets_;

    uint64_t smallest_rotation(uint64_t support) const noexcept
    {
        uint64_t smallest = support;
        for (size_t i = 1; i < dimension_; ++i) {
            support = rotate_right(support, dimension_);
            smallest = std::min(smallest, support);
        }
        return smallest;
    }

    void emit_necklace()
    {
        uint64_t support = 0;
        for (size_t i = 1; i <= dimension_; ++i) {
            if (word_[i] != 0) {
                support |= 1ULL << (dimension_ - i);
            }
        }

        const uint64_t mirrored = smallest_rotation(reflect(support, dimension_));
        if (support <= mirrored) {
            bracelets_.push_back(support);
        }
    }

    void generate_necklaces(size_t position, size_t period, size_t ones)
    {
        if (position > dimension_) {
            if (dimension_ % period == 0 && ones == support_size_) {
                emit_necklace();
            }
            return;
        }
        if (ones > support_size_
            || ones + dimension_ - position + 1 < support_size_) {
            return;
        }

        word_[position] = word_[position - period];
        generate_necklaces(position + 1, period, ones + word_[position]);

        if (word_[position - period] == 0) {
            word_[position] = 1;
            generate_necklaces(position + 1, position, ones + 1);
        }
    }
};

// Karim-Sawada-Alamgir-Husnine BraceletFC specialized to two symbols. Unlike
// the FKM path, reversal order is tracked during recursion, so rejected mirror
// necklaces are never completed and rescanned.
class DirectFixedContentBraceletGenerator {
public:
    explicit DirectFixedContentBraceletGenerator(size_t dimension)
        : dimension_(static_cast<int>(dimension))
    {}

    const std::vector<uint64_t>& generate(size_t support_size)
    {
        bracelets_.clear();
        if (support_size == 0 || support_size == static_cast<size_t>(dimension_)) {
            bracelets_.push_back(support_size == 0
                    ? 0
                    : (1ULL << dimension_) - 1ULL);
            return bracelets_;
        }

        reset(static_cast<int>(support_size));
        word_[1] = 1;
        --remaining_[1];
        if (remaining_[1] == 0) {
            list_remove(1);
        }
        blocks_[0].symbol = 0;
        update_run_length(1);
        generate_bracelets(2, 1, 1, 2, 1, false);
        return bracelets_;
    }

private:
    struct Link {
        int next = 0;
        int previous = 0;
    };

    struct Block {
        int symbol = 0;
        int length = 0;
    };

    int dimension_;
    int support_symbol_ = 0;
    int head_ = 0;
    int block_count_ = 0;
    std::array<Link, 4> available_;
    std::array<Block, kMaxDimension + 2> blocks_;
    std::array<int, 3> remaining_;
    std::array<int, kMaxDimension + 2> word_;
    std::array<int, kMaxDimension + 2> run_;
    std::vector<uint64_t> bracelets_;

    void reset(int support_size)
    {
        std::fill(remaining_.begin(), remaining_.end(), 0);
        std::fill_n(run_.begin(), dimension_ + 2, 0);
        std::fill_n(word_.begin(), dimension_ + 1, 2);
        blocks_[0] = {};

        const int non_support_size = dimension_ - support_size;
        // The paper's CAT bound requires the most frequent symbol to be last.
        // Relabeling symbols changes only the representative, not its orbit.
        if (support_size <= non_support_size) {
            remaining_[1] = support_size;
            remaining_[2] = non_support_size;
            support_symbol_ = 1;
        } else {
            remaining_[1] = non_support_size;
            remaining_[2] = support_size;
            support_symbol_ = 2;
        }

        for (int symbol = 3; symbol >= 0; --symbol) {
            available_[symbol].next = symbol - 1;
            available_[symbol].previous = symbol + 1;
        }
        head_ = 2;
        block_count_ = 0;
    }

    void list_remove(int symbol) noexcept
    {
        if (symbol == head_) {
            head_ = available_[symbol].next;
        }
        const int previous = available_[symbol].previous;
        const int next = available_[symbol].next;
        available_[previous].next = next;
        available_[next].previous = previous;
    }

    void list_add(int symbol) noexcept
    {
        const int previous = available_[symbol].previous;
        const int next = available_[symbol].next;
        available_[next].previous = symbol;
        available_[previous].next = symbol;
        if (previous == 3) {
            head_ = symbol;
        }
    }

    void update_run_length(int symbol) noexcept
    {
        if (blocks_[block_count_].symbol == symbol) {
            ++blocks_[block_count_].length;
        } else {
            ++block_count_;
            blocks_[block_count_] = { symbol, 1 };
        }
    }

    void restore_run_length() noexcept
    {
        if (blocks_[block_count_].length == 1) {
            --block_count_;
        } else {
            --blocks_[block_count_].length;
        }
    }

    // Compare the run-length encoding with its reversal. A negative result
    // proves this prefix cannot become the smallest member of a bracelet.
    int check_reverse() const noexcept
    {
        int block = 1;
        while (block <= block_count_ / 2
               && blocks_[block].length == blocks_[block_count_ - block + 1].length
               && blocks_[block].symbol == blocks_[block_count_ - block + 1].symbol) {
            ++block;
        }
        if (block > block_count_ / 2) {
            return 0;
        }
        if (blocks_[block].symbol < blocks_[block_count_ - block + 1].symbol) {
            return 1;
        }
        if (blocks_[block].symbol > blocks_[block_count_ - block + 1].symbol) {
            return -1;
        }
        if (blocks_[block].length < blocks_[block_count_ - block + 1].length
            && blocks_[block + 1].symbol < blocks_[block_count_ - block + 1].symbol) {
            return 1;
        }
        if (blocks_[block].length > blocks_[block_count_ - block + 1].length
            && blocks_[block].symbol < blocks_[block_count_ - block].symbol) {
            return 1;
        }
        return -1;
    }

    void emit(int period)
    {
        if (dimension_ % period != 0) {
            return;
        }

        uint64_t support = 0;
        for (int position = 1; position <= dimension_; ++position) {
            if (word_[position] == support_symbol_) {
                support |= 1ULL << (dimension_ - position);
            }
        }
        bracelets_.push_back(support);
    }

    void generate_bracelets(
        int position,
        int period,
        int palindrome_prefix,
        int high_run_start,
        int prefix_blocks,
        bool reverse_smaller)
    {
        // Incrementally compare the suffix after the palindromic prefix with
        // its reversal; this replaces a full O(n) comparison at every leaf.
        if (position - 1 > (dimension_ - palindrome_prefix) / 2 + palindrome_prefix) {
            const int mirrored_position = dimension_ - position + 2 + palindrome_prefix;
            if (word_[position - 1] > word_[mirrored_position]) {
                reverse_smaller = false;
            } else if (word_[position - 1] < word_[mirrored_position]) {
                reverse_smaller = true;
            }
        }

        const int positions_left = dimension_ - position + 1;
        if (remaining_[2] == positions_left) {
            if (remaining_[2] > run_[position - period]) {
                period = dimension_;
            }
            if (remaining_[2] > 0 && position != palindrome_prefix + 1
                && blocks_[prefix_blocks + 1].symbol == 2
                && blocks_[prefix_blocks + 1].length > remaining_[2]) {
                reverse_smaller = true;
            }
            if (remaining_[2] > 0 && position != palindrome_prefix + 1
                && (blocks_[prefix_blocks + 1].symbol != 2
                    || blocks_[prefix_blocks + 1].length < remaining_[2])) {
                reverse_smaller = false;
            }
            if (!reverse_smaller) {
                emit(period);
            }
            return;
        }

        // Appending only low symbols cannot complete a necklace.
        if (remaining_[1] == positions_left) {
            return;
        }

        int symbol = head_;
        while (symbol >= word_[position - period]) {
            run_[high_run_start] = position - high_run_start;
            update_run_length(symbol);
            --remaining_[symbol];
            if (remaining_[symbol] == 0) {
                list_remove(symbol);
            }
            word_[position] = symbol;

            const int next_high_run_start = symbol == 2 ? high_run_start : position + 1;
            const int next_period = symbol == word_[position - period] ? period : position;
            const int reverse_order = check_reverse();
            if (reverse_order == 0) {
                generate_bracelets(
                    position + 1, next_period, position, next_high_run_start,
                    block_count_, false);
            } else if (reverse_order == 1) {
                generate_bracelets(
                    position + 1, next_period, palindrome_prefix,
                    next_high_run_start, prefix_blocks, reverse_smaller);
            }

            if (remaining_[symbol] == 0) {
                list_add(symbol);
            }
            ++remaining_[symbol];
            restore_run_length();
            symbol = available_[symbol].next;
        }
        word_[position] = 2;
    }
};

template <typename Generator>
std::pair<size_t, uint64_t> generate_dimension(size_t dimension)
{
    Generator generator(dimension);
    size_t count = 0;
    uint64_t checksum = 0;
    for (size_t support_size = 1; support_size <= dimension; ++support_size) {
        const auto& bracelets = generator.generate(support_size);
        count += bracelets.size();
        if (!bracelets.empty()) {
            checksum = checksum * 1315423911ULL + bracelets.front();
            checksum = checksum * 1315423911ULL + bracelets.back();
        }
    }
    return { count, checksum };
}

bool verify(size_t maximum_dimension)
{
    for (size_t dimension = 1; dimension <= maximum_dimension; ++dimension) {
        FkmBraceletGenerator fkm(dimension);
        DirectFixedContentBraceletGenerator direct(dimension);
        size_t total = 0;

        for (size_t support_size = 1; support_size <= dimension; ++support_size) {
            std::vector<uint64_t> expected = fkm.generate(support_size);
            std::vector<uint64_t> actual = direct.generate(support_size);
            for (uint64_t& support : expected) {
                support = canonical_bracelet(support, dimension);
            }
            for (uint64_t& support : actual) {
                support = canonical_bracelet(support, dimension);
            }
            std::sort(expected.begin(), expected.end());
            std::sort(actual.begin(), actual.end());

            if (std::adjacent_find(actual.begin(), actual.end()) != actual.end()
                || actual != expected) {
                std::cerr << "Mismatch at n=" << dimension
                          << " k=" << support_size
                          << " fkm=" << expected.size()
                          << " direct=" << actual.size() << '\n';
                return false;
            }
            total += actual.size();
        }
        std::cout << "verified n=" << dimension << " bracelets=" << total << '\n';
    }
    return true;
}

volatile uint64_t benchmark_sink = 0;

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

template <typename Generator>
int benchmark(size_t dimension, double target_seconds)
{
    auto execute_once = [dimension]() {
        const auto [count, checksum] = generate_dimension<Generator>(dimension);
        benchmark_sink = benchmark_sink * 1315423911ULL + checksum + count;
        return count;
    };
    auto measure = [&](size_t repetitions) {
        const auto start = Clock::now();
        for (size_t i = 0; i < repetitions; ++i) {
            execute_once();
        }
        return std::chrono::duration<double, std::nano>(Clock::now() - start).count();
    };

    const double target_ns = target_seconds * 1.0e9;
    const double pilot_ns = std::max(1.0, measure(1));
    const double calibration_target_ns = std::min(20.0e6, target_ns * 0.05);
    const size_t calibration_repetitions = static_cast<size_t>(
        std::max(1.0, std::ceil(calibration_target_ns / pilot_ns)));
    const double calibrated_ns
        = measure(calibration_repetitions) / static_cast<double>(calibration_repetitions);
    const size_t measured_repetitions = static_cast<size_t>(
        std::max(1.0, std::ceil(target_ns / calibrated_ns)));
    const size_t batch_count = std::min<size_t>(7, measured_repetitions);
    const size_t base = measured_repetitions / batch_count;
    const size_t remainder = measured_repetitions % batch_count;

    std::vector<double> samples;
    samples.reserve(batch_count);
    for (size_t batch = 0; batch < batch_count; ++batch) {
        const size_t repetitions = base + (batch < remainder ? 1 : 0);
        samples.push_back(measure(repetitions) / static_cast<double>(repetitions));
    }

    const auto [count, checksum] = generate_dimension<Generator>(dimension);
    std::cout << std::setprecision(17)
              << "dimension=" << dimension << '\n'
              << "bracelets=" << count << '\n'
              << "repetitions=" << measured_repetitions << '\n'
              << "batches=" << batch_count << '\n'
              << "pilot_ns=" << pilot_ns << '\n'
              << "calibrated_ns=" << calibrated_ns << '\n'
              << "median_ns=" << median(samples) << '\n'
              << "checksum=" << checksum << '\n';
    return 0;
}

size_t parse_dimension(std::string_view argument)
{
    size_t dimension = 0;
    const auto parsed = std::from_chars(
        argument.data(), argument.data() + argument.size(), dimension);
    if (parsed.ec != std::errc {} || parsed.ptr != argument.data() + argument.size()
        || dimension == 0 || dimension > kMaxDimension) {
        throw std::invalid_argument("dimension must be from 1 through 63");
    }
    return dimension;
}

void usage(const char* program)
{
    std::cerr << "Usage:\n"
              << "  " << program << " verify [MAX_DIMENSION]\n"
              << "  " << program << " benchmark {fkm|direct} DIMENSION [TARGET_SECONDS]\n";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc >= 2 && std::string_view(argv[1]) == "verify" && argc <= 3) {
            const size_t maximum_dimension = argc == 3 ? parse_dimension(argv[2]) : 24;
            return verify(maximum_dimension) ? 0 : 1;
        }
        if (argc >= 4 && argc <= 5 && std::string_view(argv[1]) == "benchmark") {
            const size_t dimension = parse_dimension(argv[3]);
            const double target_seconds = argc == 5 ? std::stod(argv[4]) : 3.0;
            if (!(target_seconds > 0.0) || !std::isfinite(target_seconds)) {
                throw std::invalid_argument("target seconds must be finite and positive");
            }
            if (std::string_view(argv[2]) == "fkm") {
                return benchmark<FkmBraceletGenerator>(dimension, target_seconds);
            }
            if (std::string_view(argv[2]) == "direct") {
                return benchmark<DirectFixedContentBraceletGenerator>(dimension, target_seconds);
            }
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }

    usage(argv[0]);
    return 2;
}
