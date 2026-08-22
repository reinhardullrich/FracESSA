#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <fracessa/fraction.hpp>
#include <fracessa/support.hpp>
#include <fracessa/types.hpp>

#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace fracessa {

/** Result of the exact equilibrium test for one support. */
class candidate
{
public:
    explicit candidate(support::SupportContext& support_context)
        : support(support_context.make())
        , extended_support(support_context.make())
        , support_context_(&support_context)
    {}

    candidate(const candidate&) = delete;
    candidate& operator=(const candidate&) = delete;
    candidate(candidate&&) noexcept = default;
    candidate& operator=(candidate&&) noexcept = default;

    candidate clone() const
    {
        candidate result(*support_context_);
        result.copy_from(*this);
        return result;
    }

    void copy_from(const candidate& other)
    {
        candidate_id = other.candidate_id;
        vector = other.vector;
        support_context_->copy(support, other.support);
        support_size = other.support_size;
        support_context_->copy(extended_support, other.extended_support);
        extended_support_size = other.extended_support_size;
        multiplier = other.multiplier;
        is_ess = other.is_ess;
        stability = other.stability;
        payoff = other.payoff;
        payoff_dbl = other.payoff_dbl;
    }

    /// One-based stored-representative ID in deterministic search order; attempted and represented orbit members are not IDs.
    size_t candidate_id = 0;

    /// Full n-dimensional exact mixed strategy. Entries outside I(x) are zero.
    std::vector<numeric::fraction> vector;
    /// I(x), the strategies used with positive probability.
    support::Support support;
    /// Number of strategies in I(x).
    size_t support_size = 0;
    /// J(x), all pure strategies that are best replies to x.
    support::Support extended_support;
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
            if (i < vector.size() - 1) oss << ",";
        }
        oss << ";" << support_string() << ";"
            << support_size << ";"
            << extended_support_string() << ";"
            << extended_support_size << ";";
        if (multiplier) oss << *multiplier;
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
    std::string mask_to_string(const support::Support& value) const
    {
        numeric::integer result;
        result.set_string(support_context_->to_bitstring(value), 2);
        return result.to_string();
    }

    support::SupportContext* support_context_;
};

} // namespace fracessa

#endif // CANDIDATE_H
