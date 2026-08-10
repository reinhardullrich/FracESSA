#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <linalg/matrix_fraction.hpp>
#include <linalg/integer.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/bitset_multiword.hpp>
#include <optional>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <type_traits>

/**
 * Result of the exact equilibrium test for one support.
 *
 * A candidate is already a symmetric Nash equilibrium, but it is not
 * necessarily evolutionarily stable. `support` is I(x), the strategies used
 * with positive probability. `extended_support` is J(x), all pure strategies
 * that are best replies to x; it can therefore be larger than I(x).
 * `check_stability()` supplies the final ESS decision and reason code.
 */
template<class SupportMask>
class basic_candidate
{
    public:
        basic_candidate() = default;
        explicit basic_candidate(size_t dimension)
            : support(make_empty_support(dimension)), extended_support(make_empty_support(dimension))
        {}

        // IDs count stored representatives in deterministic search order, not
        // every support that was attempted or represented circular variant.
        size_t candidate_id = 0;

        // Full n-dimensional mixed strategy. Entries outside I(x) are zero.
        linalg::matrix_frc vector;
        SupportMask support;
        size_t support_size = 0;
        SupportMask extended_support;
        size_t extended_support_size = 0;

        // Number of distinct rotations/reflections represented by this row.
        // Ordinary, non-circular candidates have no multiplier.
        std::optional<size_t> multiplier;

        bool is_ess = false;
        std::string stability;
        fraction payoff;

        // Lossy copy used only in output; all mathematical decisions use payoff.
        double payoff_dbl = 0.0;

        std::string support_string() const { return mask_to_string(support); }
        std::string extended_support_string() const { return mask_to_string(extended_support); }

        std::string to_string() const
        {
            std::ostringstream oss;
            oss << candidate_id << ";";
            for (size_t i = 0; i < vector.rows(); i++) {
                oss << vector(i, 0);
                if (i < vector.rows() - 1) {
                    oss << ",";
                }
            }
            oss << ";" << support_string() << ";"
                << support_size << ";"
                << extended_support_string() << ";"
                << extended_support_size << ";";
            if (multiplier)
                oss << *multiplier;
            oss << ";" << is_ess << ";"
                << stability << ";"
                << payoff << ";"
                << std::setprecision(std::numeric_limits<double>::max_digits10) << payoff_dbl;
            return oss.str();
        }

        static std::string header()
        {
            return "candidate_id;vector;support;support_size;extended_support;extended_support_size;multiplier;is_ess;stability;payoff;payoff_dbl";
        }

    private:
        static std::string mask_to_string(const SupportMask& value)
        {
            if constexpr (std::is_same_v<SupportMask, bitset64>) {
                return bs64::to_string(value);
            } else {
                linalg::integer result;
                result.set_string(value.to_bitstring(), 2);
                return result.to_string();
            }
        }

        static SupportMask make_empty_support(size_t dimension)
        {
            if constexpr (std::is_same_v<SupportMask, bitset64>) return 0;
            else return SupportMask(dimension);
        }
};

using candidate = basic_candidate<bitset64>;
using multiword_candidate = basic_candidate<bitset_multiword>;

#endif // CANDIDATE_H
