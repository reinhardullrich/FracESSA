#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/copositivity.hpp>
#include <rational_linalg/matrix_fraction.hpp>
#include <rational_linalg/matrix_double.hpp>
#include <rational_linalg/positive_definite_fraction.hpp>
#include <rational_linalg/positive_definite_double.hpp>
#include <iostream>

void fracessa::check_stability()
{
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

    if (extended_support_size_reduced == 0)
    {
        if (conf_with_log_)
            logger_->info("Reason: true_pure_ess");
        candidate_.stability = "T_pure_ess";
        candidate_.is_ess = true;
        return;
    }
    
    auto& Bee_double = matrix_server_.get_bee_matrix_double(extended_support_reduced, m);
    
    if (Bee_double.is_positive_definite()) {
        if (conf_with_log_)
            logger_->info("Reason: true_posdef_double");
        candidate_.stability = "T_pd_double";
        candidate_.is_ess = true;
        return;
    }
    
    auto& Bee = matrix_server_.get_bee_matrix_fraction(extended_support_reduced, m);

    if (conf_with_log_) {
        logger_->info("matrix bee:\n{}", Bee.to_log_string());
    }

    if (Bee.is_positive_definite()) {
        if (conf_with_log_)
            logger_->info("Reason: true_posdef_rational");
        candidate_.stability = "T_pd_rat";
        candidate_.is_ess = true;
        return;
    }

    bitset64 kay = bs64::subtract(candidate_.extended_support, candidate_.support);
    size_t kay_size = bs64::count_set_bits(kay);

    if (conf_with_log_)
        logger_->info("kay: {}", bs64::to_bitstring(kay, dimension_));

    if (kay_size == 0 || kay_size == 1) {
        if (conf_with_log_)
            logger_->info("Reason: false_not_posdef_and_kay_0_1");
        candidate_.stability = "F_not_pd_kay_0_1";
        candidate_.is_ess = false;
        return;
    }
   
    const bitset64 jay = extended_support_reduced;
    const bitset64 jay_minus_kay = bs64::subtract(jay, kay);
    const size_t r = bs64::count_set_bits(jay_minus_kay);

    std::vector<bitset64> kay_vee(r + 1);
    std::vector<size_t> kay_vee_size(r + 1);
    std::vector<bitset64> jay_without_kay_vee(r + 1);
    std::vector<rational_linalg::matrix_fraction> bee_vee(r + 1);

    kay_vee[0] = jay;
    kay_vee_size[0] = extended_support_size_reduced;
    jay_without_kay_vee[0] = jay_minus_kay;
    bee_vee[0] = Bee;

    if (conf_with_log_) {
        logger_->info("Partial Copositivity Check:");
        logger_->info("v=0: kay_vee={}, size={}, jay\\kay={}, r={}",
                    bs64::to_bitstring(kay_vee[0], dimension_),
                    kay_vee_size[0],
                    bs64::to_bitstring(jay_without_kay_vee[0], dimension_),
                    r);
        logger_->info("bee_vee[0]:\n{}", bee_vee[0].to_log_string());
    }

    for (size_t v = 1; v <= r; ++v) {
        const size_t iv_pos = bs64::find_pos_first_set_bit(jay_without_kay_vee[v-1]);
        
        jay_without_kay_vee[v] = bs64::subtract(jay_without_kay_vee[v-1], bs64::single_bit_at_pos(iv_pos));
        kay_vee[v] = bs64::subtract(kay_vee[v-1], bs64::single_bit_at_pos(iv_pos));
        kay_vee_size[v] = kay_vee_size[v-1] - 1;
        
        const bitset64 bits_before_iv = bs64::bits_before_pos(kay_vee[v-1], iv_pos);
        const size_t pivot_pos = bs64::count_set_bits(bits_before_iv);
        
        if (conf_with_log_) {
            logger_->info("v={}: kay_vee={}, size={}, jay\\kay={}, iv_pos={}, pivot_pos={}",
                        v,
                        bs64::to_bitstring(kay_vee[v], dimension_),
                        kay_vee_size[v],
                        bs64::to_bitstring(jay_without_kay_vee[v], dimension_),
                        static_cast<unsigned int>(iv_pos),
                        static_cast<unsigned int>(pivot_pos));
        }
        
        const fraction& pivot = bee_vee[v-1](pivot_pos, pivot_pos);
        if (pivot <= fraction::zero()) {
            if (conf_with_log_) {
                logger_->info("Reason: false_not_partial_copositive (pivot={} at pos={})", pivot.to_string(), pivot_pos);
            }
            candidate_.stability = "F_not_part_copos";
            candidate_.is_ess = false;
            return;
        }
        
        bee_vee[v] = rational_linalg::matrix_fraction(kay_vee_size[v], kay_vee_size[v]);
        const size_t n_old = kay_vee_size[v-1];
        const auto& B_old = bee_vee[v-1];
        auto& B_new = bee_vee[v];
        
        for (size_t i_old = 0, i_new = 0; i_old < n_old; ++i_old) {
            if (i_old == pivot_pos) continue;            
            for (size_t j_old = 0, j_new = 0; j_old < n_old; ++j_old) {
                if (j_old == pivot_pos) continue;               
                fraction::mul(B_new(i_new, j_new), pivot, B_old(i_old, j_old));
                B_new(i_new, j_new).submul(B_old(i_old, pivot_pos), B_old(pivot_pos, j_old));
                ++j_new;
            }
            ++i_new;
        }
        
        if (conf_with_log_) {
            logger_->info("bee_vee[{}]:\n{}", v, bee_vee[v].to_log_string());
        }
    }    

    if (conf_with_log_)
        logger_->info("Copositivity Check:");

    if (rational_linalg::isStrictlyCopositiveMemoized(bee_vee[r])) {
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
