#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <fracessa/fraction.hpp>
#include <fracessa/types.hpp>
#include <fracessa/bitset.hpp>
#include <fracessa/bitset_multiword.hpp>
#include <optional>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <type_traits>
#include <vector>

namespace fracessa {

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
        /// Construct an empty candidate for the one-word support representation.
        basic_candidate() = default;
        /// Construct an empty candidate whose support masks have the requested runtime dimension.
        explicit basic_candidate(size_t dimension)
            : support(make_empty_support(dimension)), extended_support(make_empty_support(dimension))
        {}

        /// One-based stored-representative ID in deterministic search order; attempted and represented orbit members are not IDs.
        size_t candidate_id = 0;

        /// Full n-dimensional exact mixed strategy. Entries outside I(x) are zero.
        std::vector<numeric::fraction> vector;
        /// I(x), the strategies used with positive probability.
        SupportMask support;
        /// Number of strategies in I(x).
        size_t support_size = 0;
        /// J(x), all pure strategies that are best replies to x.
        SupportMask extended_support;
        /// Number of strategies in J(x).
        size_t extended_support_size = 0;

        /// Number of distinct rotations/reflections represented by this row; ordinary candidates have no multiplier.
        std::optional<size_t> multiplier;

        /// Whether exact stability checking accepted this candidate as an ESS.
        bool is_ess = false;
        /// Machine-readable reason code produced by exact stability checking.
        std::string stability;
        /// Exact equilibrium payoff x^T A x.
        numeric::fraction payoff;

        /// Lossy output-only approximation of `payoff`; no mathematical decision uses it.
        double payoff_dbl = 0.0;

        /// Return `support` as an arbitrary-width decimal integer bit mask.
        std::string support_string() const { return mask_to_string(support); }
        /// Return `extended_support` as an arbitrary-width decimal integer bit mask.
        std::string extended_support_string() const { return mask_to_string(extended_support); }

        /// Serialize this candidate as one semicolon-separated CLI row.
        std::string to_string() const
        {
            std::ostringstream oss;
            oss << candidate_id << ";";
            for (size_t i = 0; i < vector.size(); i++) {
                oss << vector[i];
                if (i < vector.size() - 1) {
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

        /// Return the semicolon-separated column names produced by `to_string()`.
        static std::string header()
        {
            return "candidate_id;vector;support;support_size;extended_support;extended_support_size;multiplier;is_ess;stability;payoff;payoff_dbl";
        }

    private:
        static std::string mask_to_string(const SupportMask& value)
        {
            if constexpr (std::is_same_v<SupportMask, support::bitset>) {
                return support::to_string(value);
            } else {
                numeric::integer result;
                result.set_string(value.to_bitstring(), 2);
                return result.to_string();
            }
        }

        static SupportMask make_empty_support(size_t dimension)
        {
            if constexpr (std::is_same_v<SupportMask, support::bitset>) return 0;
            else return SupportMask(dimension);
        }
};

using candidate = basic_candidate<support::bitset>;
using candidate_multiword = basic_candidate<support::bitset_multiword>;

} // namespace fracessa

#endif // CANDIDATE_H
