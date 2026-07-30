#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <linalg/matrix_fraction.hpp>
#include <fracessa/bitset64.hpp>
#include <optional>
#include <string>
#include <sstream>
#include <iomanip>

/*
 * Result of the exact equilibrium test for one support.
 *
 * A candidate is already a symmetric Nash equilibrium, but it is not
 * necessarily evolutionarily stable. `support` is I(x), the strategies used
 * with positive probability. `extended_support` is J(x), all pure strategies
 * that are best replies to x; it can therefore be larger than I(x).
 * `check_stability()` supplies the final ESS decision and reason code.
 */
class candidate
{
    public:
        // IDs count stored representatives in deterministic search order, not
        // every support that was attempted or represented circular variant.
        size_t candidate_id = 0;

        // Full n-dimensional mixed strategy. Entries outside I(x) are zero.
        linalg::matrix_frc vector;
        bitset64 support;
        size_t support_size = 0;
        bitset64 extended_support;
        size_t extended_support_size = 0;

        // Number of distinct rotations/reflections represented by this row.
        // Ordinary, non-circular candidates have no multiplier.
        std::optional<size_t> multiplier;

        bool is_ess = false;
        std::string stability;
        fraction payoff;

        // Lossy copy used only in output; all mathematical decisions use payoff.
        double payoff_dbl = 0.0;

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
            oss << ";" << bs64::to_string(support) << ";"
                << support_size << ";"
                << bs64::to_string(extended_support) << ";"
                << extended_support_size << ";";
            if (multiplier)
                oss << *multiplier;
            oss << ";" << is_ess << ";"
                << stability << ";"
                << payoff << ";"
                << std::fixed << std::setprecision(6) << payoff_dbl;
            return oss.str();
        }

        static std::string header()
        {
            return "candidate_id;vector;support;support_size;extended_support;extended_support_size;multiplier;is_ess;stability;payoff;payoff_dbl";
        }
};

#endif // CANDIDATE_H
