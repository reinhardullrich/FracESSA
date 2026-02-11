#pragma once

#include <flint/flint.h>  // For slong and basic types
#include <flint/fmpq.h>
#include <stdexcept>
#include <string>
#include <ostream>
#include <cstdlib>  // For free()

namespace linalg {

/*
 * Thin RAII wrapper around FLINT's fmpq_t rational type.
 *
 * Why this exists:
 * - Keep exact arithmetic semantics (no rounding) in ESS verification stages.
 * - Expose hot operations (addmul/submul/div) in-place to limit temporary objects.
 * - Provide C++ value semantics while preserving FLINT canonicalization behavior.
 */
class fraction {
private:
    fmpq_t data_;
    
public:
    // Default constructor (zero)
    fraction() noexcept {
        fmpq_init(data_);
    }
    
    // From long
    explicit fraction(long num, long den = 1) noexcept {
        fmpq_init(data_);
        fmpq_set_si(data_, num, den);
        fmpq_canonicalise(data_);
    }
    
    // From long long
    explicit fraction(long long num, long long den = 1) noexcept {
        fmpq_init(data_);
        fmpq_set_si(data_, static_cast<slong>(num), static_cast<slong>(den));
        fmpq_canonicalise(data_);
    }
    
    // From int (non-explicit so templated code can write T(0), T(1) idioms).
    fraction(int num, int den = 1) noexcept {
        fmpq_init(data_);
        fmpq_set_si(data_, static_cast<slong>(num), static_cast<slong>(den));
        fmpq_canonicalise(data_);
    }
    
    // Copy constructor
    fraction(const fraction& other) noexcept {
        fmpq_init(data_);
        fmpq_set(data_, other.data_);
    }
    
    // Move transfers ownership of the FLINT payload in O(1).
    fraction(fraction&& other) noexcept {
        fmpq_init(data_);
        fmpq_swap(data_, other.data_);
    }
    
    // Destructor
    ~fraction() noexcept {
        fmpq_clear(data_);
    }
    
    // Copy assignment
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
    
    // Direct FLINT access for low-level kernels (no wrapper overhead).
    fmpq_t& data() noexcept { return data_; }
    const fmpq_t& data() const noexcept { return data_; }
    
    // In-place arithmetic used in elimination/factorization loops.
    void add_inplace(const fraction& other) noexcept {
        fmpq_add(data_, data_, other.data_);
    }
    
    void sub_inplace(const fraction& other) noexcept {
        fmpq_sub(data_, data_, other.data_);
    }
    
    void mul_inplace(const fraction& other) noexcept {
        fmpq_mul(data_, data_, other.data_);
    }
    
    void div_inplace(const fraction& other) {
        if (fmpq_is_zero(other.data_)) {
            throw std::domain_error("Division by zero");
        }
        fmpq_div(data_, data_, other.data_);
    }
    
    // Combined operations used heavily in matrix kernels: data_ += a*b, data_ -= a*b.
    void addmul(const fraction& a, const fraction& b) noexcept {
        fmpq_addmul(data_, a.data_, b.data_);
    }
    
    void submul(const fraction& a, const fraction& b) noexcept {
        fmpq_submul(data_, a.data_, b.data_);
    }

    // Destination-first helpers avoid returning temporary fraction objects.
    static void mul(fraction& res, const fraction& a, const fraction& b) noexcept {
        fmpq_mul(res.data_, a.data_, b.data_);
    }

    static void div(fraction& res, const fraction& a, const fraction& b) {
        if (fmpq_is_zero(b.data_)) {
            throw std::domain_error("Division by zero");
        }
        fmpq_div(res.data_, a.data_, b.data_);
    }

    static void add(fraction& res, const fraction& a, const fraction& b) noexcept {
        fmpq_add(res.data_, a.data_, b.data_);
    }

    static void sub(fraction& res, const fraction& a, const fraction& b) noexcept {
        fmpq_sub(res.data_, a.data_, b.data_);
    }

    int sgn() const noexcept {
        return fmpq_sgn(data_);
    }
    
    // Arithmetic operators stay available for readability in non-critical code paths.
    
    fraction operator-(const fraction& other) const {
        fraction result;
        fmpq_sub(result.data_, data_, other.data_);
        return result;
    }
    
    fraction operator*(const fraction& other) const {
        fraction result;
        fmpq_mul(result.data_, data_, other.data_);
        return result;
    }
    
    fraction operator/(const fraction& other) const {
        if (fmpq_is_zero(other.data_)) {
            throw std::domain_error("Division by zero");
        }
        fraction result;
        fmpq_div(result.data_, data_, other.data_);
        return result;
    }
    
    fraction operator-() const {
        fraction result;
        fmpq_neg(result.data_, data_);
        return result;
    }
    
    // Compound assignment forwards to in-place FLINT operations.
    fraction& operator+=(const fraction& other) noexcept {
        add_inplace(other);
        return *this;
    }
    
    fraction& operator-=(const fraction& other) noexcept {
        sub_inplace(other);
        return *this;
    }
    
    fraction& operator*=(const fraction& other) noexcept {
        mul_inplace(other);
        return *this;
    }
    
    fraction& operator/=(const fraction& other) {
        div_inplace(other);
        return *this;
    }
    
    // Exact comparisons on canonicalized rationals.
    bool operator==(const fraction& other) const noexcept {
        return fmpq_equal(data_, other.data_);
    }
    
    bool operator!=(const fraction& other) const noexcept {
        return !fmpq_equal(data_, other.data_);
    }
    
    bool operator<(const fraction& other) const noexcept {
        return fmpq_cmp(data_, other.data_) < 0;
    }
    
    bool operator<=(const fraction& other) const noexcept {
        return fmpq_cmp(data_, other.data_) <= 0;
    }
    
    bool operator>(const fraction& other) const noexcept {
        return fmpq_cmp(data_, other.data_) > 0;
    }
    
    bool operator>=(const fraction& other) const noexcept {
        return fmpq_cmp(data_, other.data_) >= 0;
    }
    
    bool is_zero() const noexcept {
        return fmpq_is_zero(data_);
    }

    // Lossy conversion, used only for reporting and fast prefilters.
    double to_dbl() const noexcept {
        return fmpq_get_d(data_);
    }
    
    // String representation
    std::string to_string() const {
        char* str = fmpq_get_str(nullptr, 10, data_);
        if (str == nullptr) {
            return "0";
        }
        std::string result(str);
        free(str);  // FLINT's fmpq_get_str uses malloc, so use free()
        return result;
    }
    
    // Stream output writes FLINT string directly and frees FLINT-allocated memory.
    friend std::ostream& operator<<(std::ostream& os, const fraction& r) {
        char* str = fmpq_get_str(nullptr, 10, r.data_);
        if (str == nullptr) {
            os << "0";
        } else {
            os << str;  // Write directly to stream, no string allocation
            free(str);
        }
        return os;
    }
    // Shared immutable constants avoid repeated tiny object construction.
    static const fraction& zero() noexcept {
        static const fraction z(0);
        return z;
    }
    
    static const fraction& one() noexcept {
        static const fraction o(1);
        return o;
    }
    
    static const fraction& neg_one() noexcept {
        static const fraction n(-1);
        return n;
    }
    
    static const fraction& two() noexcept {
        static const fraction t(2);
        return t;
    }
};

} // namespace linalg

// Global type alias for convenience
typedef linalg::fraction fraction;
