#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/copositive_fraction.hpp>
#include <linalg/matrix_fraction.hpp>
#include <stdexcept>
#include <utility>

/*
 * Exact stability test for a candidate equilibrium p.
 *
 * Write I=I(p) for its support and J=J(p) for all pure best replies to p.
 * Choose one reference strategy m in I and form Bomze's matrix B on J without
 * m. Let K=J\I be the unused best replies. Bomze (1992), Theorems 3.2 and 3.3,
 * say that p is an ESS exactly when y^T B y is positive for every nonzero y
 * whose coordinates in K are nonnegative; coordinates belonging to I without
 * m are unrestricted.
 *
 * The code takes the cheapest exact route that settles this condition:
 * 1) J without m is empty: p is a pure strict equilibrium;
 * 2) B is positive definite: it is positive on every nonzero vector;
 * 3) eliminate the unrestricted coordinates with Bomze's rank-one recurrence;
 * 4) test ordinary strict copositivity on the remaining K-by-K matrix.
 */

void fracessa::check_stability()
{
    // Any m in I is valid. The lowest set bit gives a deterministic choice and
    // lets the bit mask and its original strategy index be obtained cheaply.
    bitset64 bitsetm = bs64::lowest_set_bit_as_bit(candidate_.support);
    bitset64 extended_support_reduced = bs64::subtract(candidate_.extended_support, bitsetm);
    size_t m = bs64::find_pos_first_set_bit(candidate_.support);
    size_t extended_support_size_reduced = candidate_.extended_support_size - 1;

    if (conf_with_log_) {
        logger_->info("Support: {}", bs64::to_bitstring(candidate_.support, dimension_));
        logger_->info("Support size: {}", bs64::count_set_bits(candidate_.support));
        logger_->info("Extended support: {}", bs64::to_bitstring(candidate_.extended_support, dimension_));
        logger_->info("Extended support size: {}", candidate_.extended_support_size);
        logger_->info("Extended support reduced: {}", bs64::to_bitstring(extended_support_reduced, dimension_));
        logger_->info("index m: {}", m);
    }

    // J={m}: no other pure strategy ties p's payoff, so the pure equilibrium is strict.
    if (extended_support_size_reduced == 0)
    {
        if (conf_with_log_)
            logger_->info("Reason: true_pure_ess");
        candidate_.stability = "T_pure_ess";
        candidate_.is_ess = true;
        return;
    }
    
    // Positive definiteness is stronger than the required cone condition and is
    // therefore a sufficient shortcut. Rational LDL^T makes it an exact certificate.
    uint8_t reduced_indices[bs64::kMaxBitsetDimension];
    const size_t bee_size = bs64::extract_set_indices(
        extended_support_reduced, dimension_, reduced_indices);
    if (bee_matrix_.rows() != bee_size) {
        bee_matrix_ = linalg::matrix_frc(bee_size, bee_size);
    }

    const fraction const_term = fraction::two() * game_matrix_(m, m);
    for (size_t row = 0; row < bee_size; ++row) {
        const size_t i = reduced_indices[row];
        for (size_t column = 0; column <= row; ++column) {
            const size_t j = reduced_indices[column];
            fraction& value = bee_matrix_(row, column);
            fraction::add(value, game_matrix_(m, j), game_matrix_(j, m));
            value += game_matrix_(i, m);
            value += game_matrix_(m, i);
            value -= game_matrix_(i, j);
            value -= game_matrix_(j, i);
            value -= const_term;

            if (row != column) {
                bee_matrix_(column, row) = value;
            }
        }
    }
    auto& Bee = bee_matrix_;

    if (conf_with_log_) {
        logger_->info("matrix bee:\n{}", Bee.to_log_string());
    }

    if (Bee.is_positive_definite()) {
        if (conf_with_log_)
            logger_->info("Reason: true_posdef_frc");
        candidate_.stability = "T_pd_frc";
        candidate_.is_ess = true;
        return;
    }

    // K=J\I contains best replies that receive zero probability in p.
    bitset64 kay = bs64::subtract(candidate_.extended_support, candidate_.support);
    size_t kay_size = bs64::count_set_bits(kay);

    if (conf_with_log_)
        logger_->info("kay: {}", bs64::to_bitstring(kay, dimension_));

    /*
     * If K is empty, the cone is all of R^(J\{m}); strict positivity is exactly
     * positive definiteness, which already failed. The same is true for one
     * constrained coordinate: for every vector, either it or its negative has
     * that coordinate nonnegative, and both have the same quadratic value.
     */
    if (kay_size == 0 || kay_size == 1) {
        if (conf_with_log_)
            logger_->info("Reason: false_not_posdef_and_kay_0_1");
        candidate_.stability = "F_not_pd_kay_0_1";
        candidate_.is_ess = false;
        return;
    }
   
    // The paper calls the reduced extended support J. J\K is exactly I\{m},
    // the set of unrestricted coordinates that the recurrence must eliminate.
    const bitset64 jay = extended_support_reduced;
    const bitset64 jay_minus_kay = bs64::subtract(jay, kay);
    const size_t r = bs64::count_set_bits(jay_minus_kay);
    if (r != candidate_.support_size - 1) {
        throw std::runtime_error("Invariant violation: r != candidate_.support_size - 1");
    }

    // State v=0 is the original reduced problem. The current index set initially
    // contains all of J, and the remaining mask contains I\{m}, the unrestricted
    // coordinates still to remove. Only the previous reduction is needed to form
    // the next one, so two rolling matrices replace the full reduction history.
    bitset64 current_indices = jay;
    bitset64 remaining_unrestricted = jay_minus_kay;
    size_t current_size = extended_support_size_reduced;
    linalg::matrix_frc previous_matrix = Bee;

    if (conf_with_log_) {
        logger_->info("Partial Copositivity Check:");
        logger_->info("v=0: kay_vee={}, size={}, jay\\kay={}, r={}", bs64::to_bitstring(current_indices, dimension_), current_size, bs64::to_bitstring(remaining_unrestricted, dimension_), r);
        logger_->info("bee_vee[0]:\n{}", previous_matrix.to_log_string());
    }

    // Eliminate the r unrestricted coordinates one at a time. What remains is
    // indexed only by K, where the ordinary nonnegative copositivity test applies.
    for (size_t v = 1; v <= r; ++v) {
        const size_t iv_pos = bs64::find_pos_first_set_bit(remaining_unrestricted);

        // iv_pos is an original strategy index. The matrix is compact, so count
        // surviving lower indices to find the corresponding row and column.
        const bitset64 bits_before_iv = bs64::bits_before_pos(current_indices, iv_pos);
        const size_t pivot_pos = bs64::count_set_bits(bits_before_iv);

        const size_t previous_size = current_size;
        remaining_unrestricted = bs64::subtract(remaining_unrestricted, bs64::single_bit_at_pos(iv_pos));
        current_indices = bs64::subtract(current_indices, bs64::single_bit_at_pos(iv_pos));
        --current_size;

        if (conf_with_log_) {
            logger_->info("v={}: kay_vee={}, size={}, jay\\kay={}, iv_pos={}, pivot_pos={}", v, bs64::to_bitstring(current_indices, dimension_), current_size, bs64::to_bitstring(remaining_unrestricted, dimension_), static_cast<unsigned int>(iv_pos), static_cast<unsigned int>(pivot_pos));
        }

        // A positive diagonal pivot is condition (a) of Bomze's recurrence. A
        // nonpositive one proves that the required cone positivity has failed.
        const fraction& pivot = previous_matrix(pivot_pos, pivot_pos);
        if (pivot <= fraction::zero()) {
            if (conf_with_log_) {
                logger_->info("Reason: false_not_partial_copositive (pivot={} at pos={})", pivot.to_string(), pivot_pos);
            }
            candidate_.stability = "F_not_part_copos";
            candidate_.is_ess = false;
            return;
        }

        /*
         * Equation (20) in Bomze (1992). A positive pivot permits removal of
         * this unrestricted coordinate. Multiplying the Schur complement by
         * the positive pivot avoids division and preserves all relevant signs:
         *
         *     B_new(i,j) = pivot*B_old(i,j)
         *                  - B_old(i,pivot)*B_old(pivot,j).
         */
        linalg::matrix_frc next_matrix(current_size, current_size);
        const auto& B_old = previous_matrix;
        auto& B_new = next_matrix;

        // Copy the compact result while omitting the eliminated row and column.
        for (size_t i_old = 0, i_new = 0; i_old < previous_size; ++i_old) {
            if (i_old == pivot_pos) continue;
            for (size_t j_old = 0, j_new = 0; j_old < previous_size; ++j_old) {
                if (j_old == pivot_pos) continue;
                fraction::mul(B_new(i_new, j_new), pivot, B_old(i_old, j_old));
                B_new(i_new, j_new).submul(B_old(i_old, pivot_pos), B_old(pivot_pos, j_old));
                ++j_new;
            }
            ++i_new;
        }

        if (conf_with_log_) {
            logger_->info("bee_vee[{}]:\n{}", v, next_matrix.to_log_string());
        }
        previous_matrix = std::move(next_matrix);
    }

    if (conf_with_log_)
        logger_->info("Copositivity Check:");

    // All unrestricted coordinates are gone; decide positivity for nonzero y>=0 on K.
    if (linalg::is_strictly_copositive(previous_matrix)) {
        if (conf_with_log_)
            logger_->info("Reason: true_copositive");
        candidate_.stability = "T_copos";
        candidate_.is_ess = true;
        return;
    } else {
        if (conf_with_log_)
            logger_->info("Reason: false_not_copositive");
        candidate_.stability = "F_not_copos";
        candidate_.is_ess = false;
        return;
    }
}
