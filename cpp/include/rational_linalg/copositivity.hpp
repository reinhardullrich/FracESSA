#ifndef RATIONAL_LINALG_COPOSITIVITY_HPP
#define RATIONAL_LINALG_COPOSITIVITY_HPP

#include <rational_linalg/matrix_fraction.hpp>
#include <rational_linalg/lu_factor_fraction.hpp>
#include <rational_linalg/adjugate_fraction.hpp>
#include <fracessa/bitset64.hpp>
#include <unordered_map>

namespace rational_linalg {

struct bitset64_hash {
    std::size_t operator()(const bitset64& bs) const noexcept {
        return bs64::hash(bs);
    }
};

class CopositivityChecker {
private:
    std::unordered_map<bitset64, bool, bitset64_hash> memo;

    bool checkRecursive(const matrix_fraction& A, const bitset64& mask) {
        auto it = memo.find(mask);
        if (it != memo.end()) return it->second;

        size_t current_dim = bs64::count_set_bits(mask);

        if (current_dim == 1) {
            size_t idx = bs64::find_pos_first_set_bit(mask);
            bool result = A(idx, idx) > fraction::zero();
            memo[mask] = result;
            return result;
        }

        bool all_good = true;
        for (size_t i = bs64::find_pos_first_set_bit(mask); i < 64; i = bs64::find_pos_next_set_bit(mask, i)) {
            bitset64 sub_mask = bs64::clear_bit_at_pos(mask, i);
            if (!checkRecursive(A, sub_mask)) {
                all_good = false;
                break;
            }
        }
        if (!all_good) {
            memo[mask] = false;
            return false;
        }

        matrix_fraction subMat = A.principal_submatrix(mask);
        lu_factor_fraction lu(subMat);
        fraction det = lu.determinant();

        if (det <= fraction::zero()) {
            matrix_fraction adj;
            if (lu.isSingular()) {
                adj = adjugate(subMat);
            } else {
                adj = lu.inverse() * det;
            }

            if (adj.all_entries_greater_zero()) {
                memo[mask] = false;
                return false;
            }
        }

        memo[mask] = true;
        return true;
    }

public:
    bool isStrictlyCopositiveMemoized(const matrix_fraction& A) {
        size_t n = A.rows();
        memo.clear();
        memo.reserve(bs64::two_to_the_power_of(n));
        bitset64 full_mask = bs64::set_all_n_bits(n);
        return checkRecursive(A, full_mask);
    }
};

inline bool isStrictlyCopositiveMemoized(const matrix_fraction& A) {
    CopositivityChecker checker;
    return checker.isStrictlyCopositiveMemoized(A);
}

} // namespace rational_linalg

#endif // RATIONAL_LINALG_COPOSITIVITY_HPP
