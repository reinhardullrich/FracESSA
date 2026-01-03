#ifndef RATIONAL_LINALG_POSITIVE_DEFINITE_FRACTION_HPP
#define RATIONAL_LINALG_POSITIVE_DEFINITE_FRACTION_HPP

#include <rational_linalg/matrix_fraction.hpp>

namespace rational_linalg {

inline bool matrix_fraction::is_positive_definite() const {
    size_t n = rows_;
    if (n == 0) return true;
    
    matrix_fraction D(n, n);
    matrix_fraction L(n, n);
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            L(i, j) = (i == j) ? fraction::one() : fraction::zero();
        }
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            fraction aSum = fraction::zero();
            for (size_t k = 0; k < j; ++k) {
                aSum += L(i, k) * L(j, k) * D(k, k);
            }
            L(i, j) = (fraction::one() / D(j, j)) * ((*this)(i, j) - aSum);
        }
        fraction bSum = fraction::zero();
        for (size_t k = 0; k < i; ++k) {
            bSum += L(i, k) * L(i, k) * D(k, k);
        }
        D(i, i) = (*this)(i, i) - bSum;
        if (D(i, i) <= fraction::zero()) {
            return false;
        }
    }
    return true;
}

inline std::vector<double> matrix_fraction::to_double_vec() const {
    std::vector<double> result;
    result.reserve(data_.size());
    for (const auto& f : data_) {
        result.push_back(f.to_double());
    }
    return result;
}

} // namespace rational_linalg

#endif // RATIONAL_LINALG_POSITIVE_DEFINITE_FRACTION_HPP
