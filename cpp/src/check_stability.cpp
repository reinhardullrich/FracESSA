#include <coposit/safe.hpp>
#include <fracessa/fracessa.hpp>
#include <fracessa/bitset.hpp>

namespace fracessa {

namespace {

inline support::bitset support_difference(support::bitset left, support::bitset right) noexcept
{
    return support::subtract(left, right);
}

inline support::bitset_multiword support_difference(support::bitset_multiword left, const support::bitset_multiword& right) noexcept
{
    left.subtract(right);
    return left;
}

} // namespace

/*
 * Exact stability test for a candidate equilibrium p.
 *
 * Write I=I(p) for its support and J=J(p) for all pure best replies. The strategies in J\I are the outside best replies. Choosing
 * one reference m in I leaves the I\{m} coordinates unrestricted and the coordinates of the outside best replies nonnegative in
 * Bomze's stability cone.
 *
 * The candidate solve has already factored the reduced Hessian H on I\{m}. Its inertia settles the common cases directly. When
 * H is negative definite and outside best replies remain, the retained factorization constructs Bomze's scaled reduced B^(r)
 * matrix through one Schur complement. Stability then reduces to ordinary strict copositivity of that smaller exact matrix.
 */

template<class SupportMask>
void basic_analyzer<SupportMask>::check_stability()
{
    // Path 1, early acceptance: J contains only the pure support strategy.
    // No other strategy ties its payoff, so this is a strict ESS.
    if (candidate_.extended_support_size == 1) {
        set_stability_result(true, "T_pure_ess");
        return;
    }

    // Outside best replies tie the candidate payoff but receive zero probability. Only their coordinates are sign constrained.
    const SupportMask outside_best_replies = support_difference(candidate_.extended_support, candidate_.support);

    // Path 2, early rejection: H is not negative definite. A support-only direction already violates strict stability,
    // so reject the candidate.
    if (!exact_candidate_solver_.reduced_hessian_is_negative_definite()) {
        set_stability_result(false, "F_reduced_hessian_not_nd");
        return;
    }

    // Path 3, early acceptance: H is negative definite and no outside best reply remains,
    // so -2H proves that the candidate is an ESS.
    if (candidate_.extended_support == candidate_.support) {
        set_stability_result(true, "T_reduced_hessian_nd");
        return;
    }
    // The remaining case has H negative definite and at least one unused tied best reply.
    exact_candidate_solver_.build_scaled_reduced_b(candidate_.support, outside_best_replies, scaled_reduced_b_);

    log_reduced_b(outside_best_replies);

    // Coposit owns the complete strict-copositivity decision, including its exact prechecks, component reduction, and final solver.
    const bool is_ess = coposit::safe::is_strictly_copositive(scaled_reduced_b_);
    set_stability_result(is_ess, is_ess ? "T_copos" : "F_not_copos");
}

template void basic_analyzer<support::bitset>::check_stability();
template void basic_analyzer<support::bitset_multiword>::check_stability();

} // namespace fracessa
