#pragma once

#include <fracessa/bitset.hpp>
#include <fracessa/bitset_multiword.hpp>

#include <cadical.hpp>

#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fracessa::support {

/** Incrementally enumerate supports of one cardinality under learned exact SAT clauses. */
class SatSupportGenerator {
public:
    explicit SatSupportGenerator(size_t dimension)
        : dimension_(dimension)
    {
        if (dimension_ == 0) throw std::invalid_argument("SAT support dimension must be positive");
        if (dimension_ > static_cast<size_t>(std::numeric_limits<int>::max() - 1))
            throw std::overflow_error("SAT variable count exceeds CaDiCaL's integer literal range");
        next_variable_ = static_cast<int>(dimension_) + 1;

        if (!solver_.configure("sat")) throw std::runtime_error("CaDiCaL lacks its satisfiable-instance configuration");
        if (!solver_.set("ilb", 2)) throw std::runtime_error("CaDiCaL lacks incremental lazy backtracking");

        std::vector<int> wires;
        size_t padded_dimension = 1;
        while (padded_dimension < dimension_) {
            if (padded_dimension > static_cast<size_t>(std::numeric_limits<int>::max()) / 2)
                throw std::overflow_error("SAT cardinality network is too large");
            padded_dimension *= 2;
        }
        wires.reserve(padded_dimension);
        for (size_t index = 0; index < dimension_; ++index) wires.push_back(variable(index));

        if (padded_dimension != dimension_) {
            const int constant_false = new_variable();
            add_clause({-constant_false});
            wires.resize(padded_dimension, constant_false);
        }

        bitonic_sort(wires, 0, wires.size(), true);
        cardinality_outputs_.assign(wires.begin(), wires.begin() + dimension_);
    }

    SatSupportGenerator(const SatSupportGenerator&) = delete;
    SatSupportGenerator& operator=(const SatSupportGenerator&) = delete;

    void start_cardinality(size_t cardinality) noexcept
    {
        assert(cardinality >= 1 && cardinality <= dimension_);
        cardinality_ = cardinality;
    }

    template<class SupportMask>
    bool take_first(SupportMask& support)
    {
        assert(cardinality_ >= 1 && cardinality_ <= dimension_);
        solver_.assume(cardinality_outputs_[cardinality_ - 1]);
        if (cardinality_ < dimension_) solver_.assume(-cardinality_outputs_[cardinality_]);

        const int status = solver_.solve();
        if (status == CaDiCaL::UNSATISFIABLE) return false;
        if (status != CaDiCaL::SATISFIABLE) throw std::runtime_error("CaDiCaL returned an inconclusive result");

        clear(support);
        size_t support_size = 0;
        for (size_t index = 0; index < dimension_; ++index) {
            if (solver_.val(variable(index)) <= 0) continue;
            set(support, index);
            ++support_size;
        }
        assert(support_size == cardinality_);
        return true;
    }

    /** Exclude every support containing `required`. */
    template<class SupportMask>
    void add_upper_cone(const SupportMask& required)
    {
        add_mask_literals(required, false);
        solver_.add(0);
    }

    /** Exclude supports containing all of `required` and none of `invaders`. */
    template<class SupportMask>
    void add_invader_interval(const SupportMask& required, const SupportMask& invaders)
    {
        for (size_t index = 0; index < dimension_; ++index) {
            const bool required_bit = contains(required, index);
            const bool invader_bit = contains(invaders, index);
            assert(!(required_bit && invader_bit));
            if (required_bit) solver_.add(-variable(index));
            else if (invader_bit) solver_.add(variable(index));
        }
        solver_.add(0);
    }

    /** Exclude only `support`; the cardinality literal makes this clause inactive in every later layer. */
    template<class SupportMask>
    void add_exact(const SupportMask& support, size_t cardinality)
    {
        assert(cardinality >= 1 && cardinality <= dimension_);
        add_mask_literals(support, false);
        if (cardinality < dimension_) solver_.add(cardinality_outputs_[cardinality]);
        solver_.add(0);
    }

private:
    static bool contains(bitset support, size_t index) noexcept
    {
        return (support & (bitset{1} << index)) != 0;
    }

    static bool contains(const bitset_multiword& support, size_t index) noexcept
    {
        return support.is_set_at_pos(index);
    }

    static void clear(bitset& support) noexcept { support = 0; }
    static void clear(bitset_multiword& support) noexcept { support.clear_all(); }
    static void set(bitset& support, size_t index) noexcept { support = set_bit_at_pos(support, index); }
    static void set(bitset_multiword& support, size_t index) noexcept { support.set_bit_at_pos(index); }

    template<class SupportMask>
    void add_mask_literals(const SupportMask& support, bool positive)
    {
        for (size_t index = 0; index < dimension_; ++index) {
            if (!contains(support, index)) continue;
            solver_.add(positive ? variable(index) : -variable(index));
        }
    }

    int variable(size_t index) const noexcept { return static_cast<int>(index) + 1; }

    int new_variable()
    {
        if (next_variable_ == std::numeric_limits<int>::max())
            throw std::overflow_error("SAT cardinality network exceeds CaDiCaL's integer literal range");
        return next_variable_++;
    }

    void add_clause(std::initializer_list<int> literals)
    {
        for (const int literal : literals) solver_.add(literal);
        solver_.add(0);
    }

    std::pair<int, int> comparator(int first, int second)
    {
        const int high = new_variable();
        const int low = new_variable();
        add_clause({-first, high});
        add_clause({-second, high});
        add_clause({first, second, -high});
        add_clause({first, -low});
        add_clause({second, -low});
        add_clause({-first, -second, low});
        return {high, low};
    }

    void compare_exchange(std::vector<int>& wires, size_t first, size_t second, bool descending)
    {
        const auto [high, low] = comparator(wires[first], wires[second]);
        wires[first] = descending ? high : low;
        wires[second] = descending ? low : high;
    }

    void bitonic_merge(std::vector<int>& wires, size_t first, size_t count, bool descending)
    {
        if (count < 2) return;
        const size_t half = count / 2;
        for (size_t index = first; index < first + half; ++index)
            compare_exchange(wires, index, index + half, descending);
        bitonic_merge(wires, first, half, descending);
        bitonic_merge(wires, first + half, half, descending);
    }

    void bitonic_sort(std::vector<int>& wires, size_t first, size_t count, bool descending)
    {
        if (count < 2) return;
        const size_t half = count / 2;
        bitonic_sort(wires, first, half, !descending);
        bitonic_sort(wires, first + half, half, descending);
        bitonic_merge(wires, first, count, descending);
    }

    size_t dimension_;
    int next_variable_ = 1;
    size_t cardinality_ = 0;
    CaDiCaL::Solver solver_;
    std::vector<int> cardinality_outputs_;
};

} // namespace fracessa::support
