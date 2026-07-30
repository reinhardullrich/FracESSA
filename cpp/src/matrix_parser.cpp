#include <fracessa/matrix_parser.hpp>

#include <charconv>
#include <cstdint>
#include <exception>
#include <iostream>
#include <vector>

namespace matrix_parser {

namespace {

constexpr size_t kMaxDimension = 63;
constexpr size_t kFastIntegerDigits = 18;

struct DecimalComponent {
    uint64_t magnitude = 0;
    size_t digits = 0;
    bool negative = false;
    bool nonzero = false;
};

bool parse_decimal_component(const char*& position, const char* end, DecimalComponent& component) noexcept
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

void append_fraction(
    std::vector<linalg::fraction>& values,
    const std::string& matrix_str,
    size_t token_start,
    size_t token_end,
    const DecimalComponent& numerator,
    const DecimalComponent& denominator)
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

bool report_value_error(const char* detail)
{
    std::cerr << "Error: Could not convert matrix values to fraction numbers!" << std::endl;
    std::cerr << "  " << detail << std::endl;
    return false;
}

} // namespace

bool parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& A, bool& is_cs)
{
    const size_t hash_pos = matrix_str.find('#');
    if (hash_pos == std::string::npos || hash_pos == 0 || hash_pos == matrix_str.length() - 1) {
        std::cerr << "Error: String for the matrix does not include '#' as a separator between dimension and matrix!" << std::endl;
        return false;
    }

    size_t n = 0;
    const char* const dimension_begin = matrix_str.data();
    const char* const dimension_end = dimension_begin + hash_pos;
    const auto [parsed_end, parse_error] = std::from_chars(dimension_begin, dimension_end, n);
    if (parse_error != std::errc{} || parsed_end != dimension_end) {
        std::cerr << "Error: The given dimension could not be converted into an integer number!" << std::endl;
        return false;
    }

    /*
     * Masks have 64 storage bits, but exhaustive search needs a representable
     * one-past-the-end value 2^n. Therefore the analyzer stops at n=63.
     */
    if (n == 0 || n > kMaxDimension) {
        std::cerr << "Error: Parser supports dimensions in [1, " << kMaxDimension << "], got " << n << std::endl;
        return false;
    }

    const size_t expected_cs = n / 2;
    const size_t expected_sym = n * (n + 1) / 2;
    std::vector<linalg::fraction> rational_values;
    rational_values.reserve(expected_sym);

    const char* const begin = matrix_str.data();
    const char* position = begin + hash_pos + 1;
    const char* const end = begin + matrix_str.size();
    try {
        while (position < end) {
            const size_t token_start = static_cast<size_t>(position - begin);
            DecimalComponent numerator;
            if (!parse_decimal_component(position, end, numerator)) {
                return report_value_error("Invalid rational numerator");
            }

            DecimalComponent denominator;
            if (position < end && *position == '/') {
                ++position;
                if (!parse_decimal_component(position, end, denominator)) {
                    return report_value_error("Invalid rational denominator");
                }
                if (!denominator.nonzero) {
                    return report_value_error("Rational denominator cannot be zero");
                }
            }

            if (position < end && *position != ',') {
                return report_value_error("Invalid character in rational value");
            }

            const size_t token_end = static_cast<size_t>(position - begin);
            append_fraction(
                rational_values, matrix_str, token_start, token_end,
                numerator, denominator);

            if (position < end) {
                ++position;
                if (position == end) {
                    return report_value_error("Empty fraction token");
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: Could not convert matrix values to fraction numbers!" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return false;
    }

    const size_t actual_size = rational_values.size();

    if (actual_size == expected_cs) {
        A = linalg::create_circular_symmetric(n, rational_values);
        is_cs = true;
    } else if (actual_size == expected_sym) {
        A = linalg::create_symmetric(n, rational_values);
        is_cs = false;
    } else {
        std::cerr << "Error: Expected " << expected_cs << " (CS) or " << expected_sym << " (Sym) values, got " << actual_size << std::endl;
        return false;
    }

    return true;
}

} // namespace matrix_parser
