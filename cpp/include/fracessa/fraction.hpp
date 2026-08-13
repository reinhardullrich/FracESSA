#pragma once

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <ostream>
#include <string>

#include <fracessa/types.hpp>

namespace fracessa::numeric {

/*
 * Owning exact rational value used only at FracESSA's public output boundary.
 *
 * Every object initializes and clears one FLINT rational. Candidate search and
 * stability remain in integer arithmetic; this type materializes accepted
 * probabilities and payoffs for comparison and serialization.
 *
 * Numeric constructors exist for concise exact-output tests. They and
 * set_ratio() require a nonzero denominator.
 */
class fraction {
private:
    fmpq_t data_;

    static void set_signed_ratio(fmpq_t target, slong num, slong den) noexcept {
        if (den < 0) {
            num = -num;
            den = -den;
        }
        fmpq_set_si(target, num, static_cast<ulong>(den));
        fmpq_canonicalise(target);
    }
    
public:
    fraction() noexcept {
        fmpq_init(data_);
    }
    
    explicit fraction(long long num, long long den = 1) noexcept {
        fmpq_init(data_);
        set_signed_ratio(data_, static_cast<slong>(num), static_cast<slong>(den));
    }

    // Implicit int construction keeps exact expected values concise in C++ tests.
    fraction(int num, int den = 1) noexcept {
        fmpq_init(data_);
        set_signed_ratio(data_, static_cast<slong>(num), static_cast<slong>(den));
    }
    
    fraction(const fraction& other) noexcept {
        fmpq_init(data_);
        fmpq_set(data_, other.data_);
    }
    
    // Move transfers ownership of the FLINT payload in O(1).
    fraction(fraction&& other) noexcept {
        fmpq_init(data_);
        fmpq_swap(data_, other.data_);
    }
    
    ~fraction() noexcept {
        fmpq_clear(data_);
    }
    
    fraction& operator=(const fraction& other) noexcept {
        if (this != &other) {
            fmpq_set(data_, other.data_);
        }
        return *this;
    }
    
    // Move-assignment mirrors move-ctor semantics with swap.
    fraction& operator=(fraction&& other) noexcept {
        if (this != &other) {
            fmpq_swap(data_, other.data_);
        }
        return *this;
    }
    
    void set_zero() noexcept {
        fmpq_zero(data_);
    }

    // Materialize an exact rational from C++ integer wrappers without exposing FLINT at the algorithm call site.
    void set_ratio(integer::const_reference numerator, integer::const_reference denominator) noexcept {
        fmpq_set_fmpz_frac(data_, numerator.native_handle(), denominator.native_handle());
    }
    
    // Exact comparisons on canonicalized rationals.
    bool operator==(const fraction& other) const noexcept {
        return fmpq_equal(data_, other.data_);
    }
    
    // Lossy conversion for reporting only, never a mathematical decision.
    double to_dbl() const noexcept {
        return fmpq_get_d(data_);
    }
    
    // FLINT allocates text buffers, so both conversion paths use flint_free.
    std::string to_string() const {
        char* str = fmpq_get_str(nullptr, 10, data_);
        if (str == nullptr) {
            return "0";
        }
        std::string result(str);
        flint_free(str);
        return result;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const fraction& r) {
        char* str = fmpq_get_str(nullptr, 10, r.data_);
        if (str == nullptr) {
            os << "0";
        } else {
            os << str;  // Avoid an additional std::string copy.
            flint_free(str);
        }
        return os;
    }
};

} // namespace fracessa::numeric
