#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/copositive_fraction.hpp>

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
    const size_t kay_size = bs64::count_set_bits(kay);

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
    find_candidate_safe_.build_scaled_reduced_b(candidate_.support, kay, scaled_reduced_b_);

    if (conf_with_log_) {
        logger_->info("kay: {}", bs64::to_bitstring(kay, dimension_));
        logger_->info("scaled reduced B matrix:\n{}", scaled_reduced_b_.to_log_string());
    }

    // Path 4, early acceptance: The scaled reduced B matrix is positive definite. This implies strict copositivity,
    // so accept the candidate as an ESS.
    if (scaled_reduced_b_.is_positive_definite()) {
        if (conf_with_log_) logger_->info("Reason: true_posdef_frc");
        candidate_.stability = "T_pd_frc";
        candidate_.is_ess = true;
        return;
    }

    // Path 5, early rejection: K has exactly one coordinate, because Path 3 excluded K=empty.
    // Failure of positive definiteness therefore also means failure of strict copositivity, so reject the candidate.
    if (kay_size <= 1) {
        if (conf_with_log_) logger_->info("Reason: false_not_posdef_and_kay_0_1");
        candidate_.stability = "F_not_pd_kay_0_1";
        candidate_.is_ess = false;
        return;
    }

    // Path 6, final decision: K has at least two coordinates. Strict copositivity of the scaled reduced B matrix
    // is now the exact remaining test; accept the candidate exactly when this test succeeds.
    if (linalg::is_strictly_copositive(scaled_reduced_b_)) {
        if (conf_with_log_) logger_->info("Reason: true_copositive");
        candidate_.stability = "T_copos";
        candidate_.is_ess = true;
    } else {
        if (conf_with_log_) logger_->info("Reason: false_not_copositive");
        candidate_.stability = "F_not_copos";
        candidate_.is_ess = false;
    }
}
