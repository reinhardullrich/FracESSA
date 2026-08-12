#pragma once

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <linalg/matrix_fraction.hpp>

namespace matrix_parser {

namespace detail {

inline constexpr size_t kFastIntegerDigits = 18;

inline size_t checked_symmetric_value_count(size_t dimension)
{
    const size_t left = dimension % 2 == 0 ? dimension / 2 : dimension;
    const size_t right = dimension % 2 == 0 ? dimension + 1 : dimension / 2 + 1;
    if (left != 0 && right > std::numeric_limits<size_t>::max() / left)
        throw std::invalid_argument("Matrix dimension is too large: " + std::to_string(dimension));
    return left * right;
}

inline void validate_dense_matrix_size(size_t dimension)
{
    const size_t max_entries = std::vector<linalg::fraction>().max_size();
    if (dimension > max_entries / dimension)
        throw std::invalid_argument("Matrix dimension is too large: " + std::to_string(dimension));
}

struct DecimalComponent {
    uint64_t magnitude = 0;
    size_t digits = 0;
    bool negative = false;
    bool nonzero = false;
};

inline bool parse_decimal_component(const char*& position, const char* end, DecimalComponent& component) noexcept
{
    if (position < end && *position == '-') {
        component.negative = true;
        ++position;
    }

    const char* const digits_begin = position;
    while (position < end && *position >= '0' && *position <= '9') {
        const unsigned digit = static_cast<unsigned>(*position - '0');
        // Unsigned overflow is defined and ignored when the text takes the arbitrary-precision path.
        component.magnitude = component.magnitude * 10 + digit;
        component.nonzero = component.nonzero || digit != 0;
        ++position;
    }
    component.digits = static_cast<size_t>(position - digits_begin);
    return component.digits != 0;
}

inline void append_fraction(std::vector<linalg::fraction>& values, const std::string& matrix_str, size_t token_start,
                            size_t token_end, const DecimalComponent& numerator, const DecimalComponent& denominator)
{
    if (numerator.digits <= kFastIntegerDigits &&
        (denominator.digits == 0 || denominator.digits <= kFastIntegerDigits)) {
        long long num = static_cast<long long>(numerator.magnitude);
        if (numerator.negative) num = -num;

        long long den = 1;
        if (denominator.digits != 0) {
            den = static_cast<long long>(denominator.magnitude);
            if (denominator.negative) den = -den;
        }
        values.emplace_back(num, den);
        return;
    }

    values.emplace_back(matrix_str.substr(token_start, token_end - token_start));
}

} // namespace detail

/**
 * Parse a validating CLI/Python `n#values` payload into an exact symmetric matrix.
 *
 * Values are comma-separated integers or exact integer fractions. Decimal notation is rejected. A list of `n*(n+1)/2` values is
 * interpreted as the upper triangle in row-major order. For `n >= 2`, a shorter list of `floor(n/2)` values is interpreted as the
 * circular-distance entries of a zero-diagonal circular-symmetric matrix.
 *
 * @param matrix_str Complete matrix payload including its positive dimension.
 * @param matrix Receives the parsed n-by-n exact rational matrix.
 * @param is_circular_symmetric Receives true for compact circular input and false for full upper-triangular input.
 * @throws std::invalid_argument if the dimension, value syntax, denominator, or number of values is invalid.
 */
inline void parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& matrix, bool& is_circular_symmetric)
{
    const size_t hash_pos = matrix_str.find('#');
    if (hash_pos == std::string::npos || hash_pos == 0 || hash_pos == matrix_str.length() - 1) {
        throw std::invalid_argument("String for the matrix does not include '#' as a separator between dimension and matrix");
    }
    if (matrix_str.find('#', hash_pos + 1) != std::string::npos) {
        throw std::invalid_argument("Multiple '#' characters found in matrix string");
    }

    size_t n = 0;
    const char* const dimension_begin = matrix_str.data();
    const char* const dimension_end = dimension_begin + hash_pos;
    const auto [parsed_end, parse_error] = std::from_chars(dimension_begin, dimension_end, n);
    if (parse_error != std::errc{} || parsed_end != dimension_end) {
        throw std::invalid_argument("The given dimension could not be converted into an integer number");
    }
    if (n == 0) throw std::invalid_argument("Matrix dimension must be positive");

    detail::validate_dense_matrix_size(n);

    const size_t expected_cs = n / 2;
    const size_t expected_sym = detail::checked_symmetric_value_count(n);
    std::vector<linalg::fraction> rational_values;
    rational_values.reserve(expected_sym < matrix_str.size() ? expected_sym : matrix_str.size());

    const char* const begin = matrix_str.data();
    const char* position = begin + hash_pos + 1;
    const char* const end = begin + matrix_str.size();
    while (position < end) {
        const size_t token_start = static_cast<size_t>(position - begin);
        detail::DecimalComponent numerator;
        if (!detail::parse_decimal_component(position, end, numerator)) {
            throw std::invalid_argument("Invalid rational numerator");
        }

        detail::DecimalComponent denominator;
        if (position < end && *position == '/') {
            ++position;
            if (!detail::parse_decimal_component(position, end, denominator)) {
                throw std::invalid_argument("Invalid rational denominator");
            }
            if (!denominator.nonzero) {
                throw std::invalid_argument("Rational denominator cannot be zero: " +
                                            matrix_str.substr(token_start,
                                                              static_cast<size_t>(position - begin) - token_start));
            }
        }

        if (position < end && *position != ',') {
            throw std::invalid_argument("Invalid character in rational value");
        }

        const size_t token_end = static_cast<size_t>(position - begin);
        detail::append_fraction(rational_values, matrix_str, token_start, token_end, numerator, denominator);

        if (position < end) {
            ++position;
            if (position == end) {
                throw std::invalid_argument("Empty fraction token");
            }
        }
    }

    const size_t actual_size = rational_values.size();
    if (actual_size == expected_cs) {
        matrix = linalg::create_circular_symmetric(n, rational_values);
        is_circular_symmetric = true;
    } else if (actual_size == expected_sym) {
        matrix = linalg::create_symmetric(n, rational_values);
        is_circular_symmetric = false;
    } else {
        throw std::invalid_argument("Expected " + std::to_string(expected_cs) + " (CS) or " +
                                    std::to_string(expected_sym) + " (Sym) values, got " + std::to_string(actual_size));
    }
}

} // namespace matrix_parser
