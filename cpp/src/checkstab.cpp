#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/copositive_integer.hpp>

/*
 * Exact stability test for a candidate equilibrium p.
 *
 * Write I=I(p) for its support, J=J(p) for all pure best replies, and K=J\I for the unused best replies. Choosing one reference
 * m in I leaves the I\{m} coordinates unrestricted and the K coordinates nonnegative in Bomze's stability cone.
 *
 * The candidate solve has already factored the reduced Hessian H on I\{m}. Its inertia settles the common cases directly. When
 * H is negative definite and K is nonempty, the retained factorization constructs Bomze's scaled reduced B^(r) matrix through one
 * Schur complement. Stability then reduces to ordinary strict copositivity of that smaller exact matrix.
 */

void fracessa::check_stability()
{
    if (conf_with_log_) {
        // Any m in I is valid. The lowest set bit matches the deterministic reference used by the exact candidate solve.
        const bitset64 bitsetm = bs64::lowest_set_bit_as_bit(candidate_.support);
        const bitset64 extended_support_reduced = bs64::subtract(candidate_.extended_support, bitsetm);
        const size_t m = bs64::find_pos_first_set_bit(candidate_.support);
        logger_->info("Support: {}", bs64::to_bitstring(candidate_.support, dimension_));
        logger_->info("Support size: {}", candidate_.support_size);
        logger_->info("Extended support: {}", bs64::to_bitstring(candidate_.extended_support, dimension_));
        logger_->info("Extended support size: {}", candidate_.extended_support_size);
        logger_->info("Extended support reduced: {}", bs64::to_bitstring(extended_support_reduced, dimension_));
        logger_->info("index m: {}", m);
    }

    // Path 1, early acceptance: J contains only the pure support strategy.
    // No other strategy ties its payoff, so this is a strict ESS.
    if (candidate_.extended_support_size == 1) {
        if (conf_with_log_) logger_->info("Reason: true_pure_ess");
        candidate_.stability = "T_pure_ess";
        candidate_.is_ess = true;
        return;
    }

    // K contains best replies that receive zero probability in p. Only their coordinates are sign constrained.
    const bitset64 kay = bs64::subtract(candidate_.extended_support, candidate_.support);

    // Path 2, early rejection: H is not negative definite. A support-only direction already violates strict stability,
    // so reject the candidate.
    if (!find_candidate_safe_.reduced_hessian_is_negative_definite()) {
        if (conf_with_log_) logger_->info("Reason: false_reduced_hessian_not_negative_definite");
        candidate_.stability = "F_reduced_hessian_not_nd";
        candidate_.is_ess = false;
        return;
    }

    // Path 3, early acceptance: H is negative definite and K is empty. No unused tied best reply remains,
    // so -2H proves that the candidate is an ESS.
    if (candidate_.extended_support == candidate_.support) {
        if (conf_with_log_) logger_->info("Reason: true_posdef_frc (from reduced Hessian)");
        candidate_.stability = "T_pd_frc";
        candidate_.is_ess = true;
        return;
    }
    // The remaining case has H negative definite and at least one unused tied best reply.
    const size_t kay_size = find_candidate_safe_.build_scaled_reduced_b(candidate_.support, kay, scaled_reduced_b_);

    if (conf_with_log_) {
        logger_->info("kay: {}", bs64::to_bitstring(kay, dimension_));
        logger_->info("scaled reduced B matrix:\n{}", scaled_reduced_b_.to_log_string());
    }

    // Path 4, direct decision: Small K has fixed exact criteria, cheaper than collecting data for the general copositivity path.
    switch (kay_size) {
        case 1:
            candidate_.is_ess = linalg::is_strictly_copositive_1x1(scaled_reduced_b_(0, 0));
            candidate_.stability = candidate_.is_ess ? "T_copos_small" : "F_not_copos_small";
            if (conf_with_log_) logger_->info("Reason: {}", candidate_.stability);
            return;
        case 2:
            candidate_.is_ess = linalg::is_strictly_copositive_2x2(
                scaled_reduced_b_(0, 0), scaled_reduced_b_(0, 1), scaled_reduced_b_(1, 1));
            candidate_.stability = candidate_.is_ess ? "T_copos_small" : "F_not_copos_small";
            if (conf_with_log_) logger_->info("Reason: {}", candidate_.stability);
            return;
        case 3:
            candidate_.is_ess = linalg::is_strictly_copositive_3x3(
                scaled_reduced_b_(0, 0), scaled_reduced_b_(0, 1), scaled_reduced_b_(0, 2),
                scaled_reduced_b_(1, 1), scaled_reduced_b_(1, 2), scaled_reduced_b_(2, 2));
            candidate_.stability = candidate_.is_ess ? "T_copos_small" : "F_not_copos_small";
            if (conf_with_log_) logger_->info("Reason: {}", candidate_.stability);
            return;
        default:
            break;
    }

    const linalg::CopositivitySignScan sign_scan = linalg::scan_copositivity_signs(scaled_reduced_b_);

    // This must be the first sign-scan decision: every other scan field is incomplete when a diagonal fails.
    // Path 5, early rejection: A nonpositive diagonal is an explicit coordinate-vector witness against strict copositivity.
    if (!sign_scan.all_diagonal_positive) {
        if (conf_with_log_) logger_->info("Reason: F_not_copos_nonpositive_diagonal");
        candidate_.stability = "F_not_copos_nonpositive_diagonal";
        candidate_.is_ess = false;
        return;
    }

    // Path 6, early acceptance: Positive diagonals and nonnegative off-diagonal entries make every term in the quadratic form
    // nonnegative, with at least one positive term for every nonzero nonnegative vector.
    if (!sign_scan.has_negative_off_diagonal) {
        if (conf_with_log_) logger_->info("Reason: T_copos_nonnegative_off_diagonal");
        candidate_.stability = "T_copos_nonnegative_off_diagonal";
        candidate_.is_ess = true;
        return;
    }

    // Path 7, early acceptance: Liqun Qi, LAA 439 (2013), Theorem 10, eq. (12): a positive diagonal plus all negative
    // off-diagonal entries in every row implies strict copositivity.
    if (sign_scan.all_negative_part_row_sums_positive) {
        if (conf_with_log_) logger_->info("Reason: T_copos_negative_part_diagonal_dominance");
        candidate_.stability = "T_copos_negative_part_diagonal_dominance";
        candidate_.is_ess = true;
        return;
    }

    // Path 8, early rejection: The all-ones vector is nonzero and nonnegative. If its exact quadratic value is nonpositive,
    // the scaled reduced B matrix is not strictly copositive.
    if (sign_scan.all_ones_quadratic_value.sign() <= 0) {
        if (conf_with_log_) logger_->info("Reason: F_not_copos_nonpositive_all_ones_value");
        candidate_.stability = "F_not_copos_nonpositive_all_ones_value";
        candidate_.is_ess = false;
        return;
    }

    // Path 9, final decision: nonnegative interactions between distinct components of the negative-entry graph cannot create a
    // nonpositive value on the nonnegative orthant. Check each component independently with the exact adaptive-cone test.
    if (linalg::CopositivityChecker::are_negative_components_strictly_copositive(
            scaled_reduced_b_, sign_scan.negative_neighbors)) {
        if (conf_with_log_) logger_->info("Reason: true_copositive");
        candidate_.stability = "T_copos";
        candidate_.is_ess = true;
    } else {
        if (conf_with_log_) logger_->info("Reason: false_not_copositive");
        candidate_.stability = "F_not_copos";
        candidate_.is_ess = false;
    }
}
