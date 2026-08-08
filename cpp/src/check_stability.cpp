#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/copositive_integer.hpp>

namespace {

inline bitset64 support_difference(bitset64 left, bitset64 right) noexcept { return bs64::subtract(left, right); }

inline bitset_multiword support_difference(bitset_multiword left, const bitset_multiword& right) noexcept
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
void basic_fracessa<SupportMask>::check_stability()
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
    const size_t outside_best_reply_count =
        exact_candidate_solver_.build_scaled_reduced_b(candidate_.support, outside_best_replies, scaled_reduced_b_);

    log_reduced_b(outside_best_replies);

    // Path 4, direct decision: One to three outside best replies have fixed exact criteria, cheaper than the general path.
    switch (outside_best_reply_count) {
        case 1:
        {
            const bool is_ess = linalg::is_strictly_copositive_1x1(scaled_reduced_b_(0, 0));
            set_stability_result(is_ess, is_ess ? "T_copos_small" : "F_not_copos_small");
            return;
        }
        case 2:
        {
            const bool is_ess = linalg::is_strictly_copositive_2x2(
                scaled_reduced_b_(0, 0), scaled_reduced_b_(0, 1), scaled_reduced_b_(1, 1));
            set_stability_result(is_ess, is_ess ? "T_copos_small" : "F_not_copos_small");
            return;
        }
        case 3:
        {
            const bool is_ess = linalg::is_strictly_copositive_3x3(
                scaled_reduced_b_(0, 0), scaled_reduced_b_(0, 1), scaled_reduced_b_(0, 2),
                scaled_reduced_b_(1, 1), scaled_reduced_b_(1, 2), scaled_reduced_b_(2, 2));
            set_stability_result(is_ess, is_ess ? "T_copos_small" : "F_not_copos_small");
            return;
        }
        default:
            break;
    }

    const auto decide_from_sign_scan = [&](const auto& sign_scan) {
        // This must be the first sign-scan decision: every other scan field is incomplete when a diagonal fails.
        // Path 5, early rejection: A nonpositive diagonal is an explicit coordinate-vector witness against strict copositivity.
        if (!sign_scan.all_diagonal_positive) {
            set_stability_result(false, "F_not_copos_nonpositive_diagonal");
            return;
        }

        // Path 6, early acceptance: Positive diagonals and nonnegative off-diagonal entries make every term in the quadratic form
        // nonnegative, with at least one positive term for every nonzero nonnegative vector.
        if (!sign_scan.has_negative_off_diagonal) {
            set_stability_result(true, "T_copos_nonnegative_off_diagonal");
            return;
        }

        // Path 7, early acceptance: Liqun Qi, LAA 439 (2013), Theorem 10, eq. (12): a positive diagonal plus all negative
        // off-diagonal entries in every row implies strict copositivity.
        if (sign_scan.all_negative_part_row_sums_positive) {
            set_stability_result(true, "T_copos_negative_part_diagonal_dominance");
            return;
        }

        // Path 8, early rejection: The all-ones vector is nonzero and nonnegative. If its exact quadratic value is nonpositive,
        // the scaled reduced B matrix is not strictly copositive.
        if (sign_scan.all_ones_quadratic_value.sign() <= 0) {
            set_stability_result(false, "F_not_copos_nonpositive_all_ones_value");
            return;
        }

        // Path 9, final decision: nonnegative interactions between distinct components of the negative-entry graph cannot create a
        // nonpositive value on the nonnegative orthant. Check each component independently with exact Hadeler enumeration.
        linalg::CopositivityChecker copositivity_checker(outside_best_reply_count);
        const bool is_ess = copositivity_checker.are_negative_components_strictly_copositive(
            scaled_reduced_b_, sign_scan.negative_neighbors);
        set_stability_result(is_ess, is_ess ? "T_copos" : "F_not_copos");
    };

    // The reduced B dimension, not the original game dimension, determines whether one or several words are needed. Keep the
    // one-word analyzer's choice at compile time so its path gains no runtime branch.
    if constexpr (std::is_same_v<SupportMask, bitset64>)
        decide_from_sign_scan(linalg::scan_copositivity_signs(scaled_reduced_b_));
    else if (outside_best_reply_count <= bs64::kMaxBitsetDimension)
        decide_from_sign_scan(linalg::scan_copositivity_signs(scaled_reduced_b_));
    else
        decide_from_sign_scan(linalg::scan_copositivity_signs_multiword(scaled_reduced_b_));
}

template void basic_fracessa<bitset64>::check_stability();
template void basic_fracessa<bitset_multiword>::check_stability();
