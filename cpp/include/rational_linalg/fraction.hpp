#pragma once

#include <flint/flint.h>  // For slong and basic types
#include <flint/fmpq.h>
#include <stdexcept>
#include <string>
#include <ostream>
#include <cstdlib>  // For free()

namespace rational_linalg {

// ============================================================================
// fraction - OPTIMIZED FLINT Rational Number Wrapper
// ============================================================================
// Thin C++17 wrapper around FLINT's fmpq_t for arbitrary-precision rationals
// ============================================================================

class fraction {
private:
    fmpq_t data_;
    
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
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
    
    // From int (non-explicit for compatibility with T(0), T(1) patterns)
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
    
    // FIXED: Correct move constructor using swap
    fraction(fraction&& other) noexcept {
        fmpq_init(data_);
        fmpq_swap(data_, other.data_);
    }
    
    // Destructor
    ~fraction() noexcept {
        fmpq_clear(data_);
    }
    
    // ========================================================================
    // Assignment Operators
    // ========================================================================
    
    // Copy assignment
    fraction& operator=(const fraction& other) noexcept {
        if (this != &other) {
            fmpq_set(data_, other.data_);
        }
        return *this;
    }
    
    // FIXED: Correct move assignment using swap
    fraction& operator=(fraction&& other) noexcept {
        if (this != &other) {
            fmpq_swap(data_, other.data_);
        }
        return *this;
    }
    
    // ========================================================================
    // Direct Access to Underlying FLINT Type (Zero Overhead)
    // ========================================================================
    
    fmpq_t& data() noexcept { return data_; }
    const fmpq_t& data() const noexcept { return data_; }
    
    // Get pointer for direct FLINT operations
    fmpq* ptr() noexcept { return data_; }
    const fmpq* ptr() const noexcept { return data_; }
    
    // ========================================================================
    // In-Place Operations (FAST - No Temporaries)
    // ========================================================================
    
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
    
    void negate_inplace() noexcept {
        fmpq_neg(data_, data_);
    }
    
    void abs_inplace() noexcept {
        fmpq_abs(data_, data_);
    }
    
    // ========================================================================
    // Arithmetic Operators (RVO-Friendly)
    // ========================================================================
    
    fraction operator+(const fraction& other) const {
        fraction result;
        fmpq_add(result.data_, data_, other.data_);
        return result;
    }
    
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
    
    // ========================================================================
    // Compound Assignment Operators (Use In-Place Operations)
    // ========================================================================
    
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
    
    // ========================================================================
    // Comparison Operators (noexcept)
    // ========================================================================
    
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
    
    // ========================================================================
    // Utility Functions
    // ========================================================================
    
    bool is_zero() const noexcept {
        return fmpq_is_zero(data_);
    }
    
    bool is_one() const noexcept {
        return fmpq_is_one(data_);
    }
    
    fraction abs() const noexcept {
        fraction result;
        fmpq_abs(result.data_, data_);
        return result;
    }
    
    fraction inverse() const {
        if (fmpq_is_zero(data_)) {
            throw std::domain_error("Cannot invert zero");
        }
        fraction result;
        fmpq_inv(result.data_, data_);
        return result;
    }
    
    // ========================================================================
    // Conversions
    // ========================================================================
    
    // Convert to double (ONLY conversion implemented)
    double to_double() const noexcept {
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
    
    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const fraction& r) {
        os << r.to_string();
        return os;
    }
};

// ============================================================================
// Helper Functions
// ============================================================================

} // namespace rational_linalg

// Global type alias for convenience
typedef rational_linalg::fraction fraction;
