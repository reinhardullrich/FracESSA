#pragma once

#include <rational_linalg/fraction.hpp>

namespace rational_linalg {

// Type traits for scalar constants
template<typename T>
struct ScalarConstants;

// Specialization for double
template<>
struct ScalarConstants<double> {
    static constexpr double zero() { return 0.0; }
    static constexpr double one() { return 1.0; }
    static constexpr double neg_one() { return -1.0; }
};

// Specialization for fraction
template<>
struct ScalarConstants<fraction> {
    static const fraction& zero() {
        static const fraction z(0);
        return z;
    }
    
    static const fraction& one() {
        static const fraction o(1);
        return o;
    }
    
    static const fraction& neg_one() {
        static const fraction n(-1);
        return n;
    }
};

// Convenience functions - return by value for double, reference for fraction
template<typename T>
inline auto zero() -> decltype(ScalarConstants<T>::zero()) {
    return ScalarConstants<T>::zero();
}

template<typename T>
inline auto one() -> decltype(ScalarConstants<T>::one()) {
    return ScalarConstants<T>::one();
}

template<typename T>
inline auto neg_one() -> decltype(ScalarConstants<T>::neg_one()) {
    return ScalarConstants<T>::neg_one();
}

} // namespace rational_linalg

