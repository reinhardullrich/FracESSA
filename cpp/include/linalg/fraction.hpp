#pragma once

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <stdexcept>
#include <string>
#include <ostream>

#include <linalg/integer.hpp>

namespace linalg {

/*
 * Owning C++ value wrapper around FLINT's fmpq_t rational type.
 *
 * Every object initializes and clears one FLINT rational. FLINT keeps values
 * canonical, so arithmetic and comparisons are exact and need no epsilon.
 * Returning operators keep formulas readable; destination-first and in-place
 * helpers avoid wrapper temporaries in the small-matrix inner loops.
 *
 * Value constructors assume every denominator is nonzero. The matrix parser
 * validates external text before construction; internal callers deliberately
 * retain this precondition to avoid repeated checks.
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

    // From FLINT-compatible base-10 rational text (`num` or `num/den`).
    explicit fraction(const std::string& value) {
        fmpq_init(data_);
        if (fmpq_set_str(data_, value.c_str(), 10) != 0) {
            fmpq_clear(data_);
            throw std::invalid_argument("Invalid rational value: " + value);
        }
        fmpq_canonicalise(data_);
    }
    
    // Implicit int construction keeps exact expressions with 0 and 1 concise.
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
    
    // Escape hatch for FLINT kernels; callers must preserve valid fmpq state.
    fmpq_t& data() noexcept { return data_; }
    const fmpq* data() const noexcept { return data_; }

    // Construct an exact rational from C++ integer wrappers without exposing FLINT at the algorithm call site.
    void set_ratio(integer::const_reference numerator, integer::const_reference denominator) noexcept {
        fmpq_set_fmpz_frac(data_, numerator.native_handle(), denominator.native_handle());
    }
    
    // In-place arithmetic used in elimination/factorization loops.
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
    
    // Returning operators favor readable formulas; hot loops use the helpers above.
    
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

    // Compound assignment writes into existing FLINT storage without a temporary.
    fraction& operator+=(const fraction& other) noexcept {
        fmpq_add(data_, data_, other.data_);
        return *this;
    }

    fraction& operator-=(const fraction& other) noexcept {
        fmpq_sub(data_, data_, other.data_);
        return *this;
    }
    
    // Used by exact LU determinant accumulation.
    fraction& operator*=(const fraction& other) noexcept {
        mul_inplace(other);
        return *this;
    }
    
    // Exact comparisons on canonicalized rationals.
    bool operator==(const fraction& other) const noexcept {
        return fmpq_equal(data_, other.data_);
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
    
    bool is_zero() const noexcept {
        return fmpq_is_zero(data_);
    }

    // Lossy conversion for reporting and fast candidate search, never an exact certificate.
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

typedef linalg::fraction fraction;
