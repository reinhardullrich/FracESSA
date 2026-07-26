#include <fracessa/matrix_parser.hpp>

#include <charconv>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace matrix_parser {

namespace {

linalg::fraction parse_fraction_token(const std::string& values, size_t start, size_t end)
{
    if (start >= end) {
        throw std::invalid_argument("empty fraction token");
    }

    const std::string token = values.substr(start, end - start);
    const size_t slash_pos = token.find('/');
    if (slash_pos != std::string::npos) {
        size_t denominator_pos = slash_pos + 1;
        if (denominator_pos < token.size() && (token[denominator_pos] == '+' || token[denominator_pos] == '-')) {
            ++denominator_pos;
        }

        if (denominator_pos == token.size()) {
            throw std::invalid_argument("Invalid rational denominator: " + token);
        }

        bool denominator_is_zero = true;
        for (size_t i = denominator_pos; i < token.size(); ++i) {
            if (token[i] < '0' || token[i] > '9') {
                throw std::invalid_argument("Invalid rational denominator: " + token);
            }
            denominator_is_zero = denominator_is_zero && token[i] == '0';
        }
        if (denominator_is_zero) {
            throw std::invalid_argument("Rational denominator cannot be zero: " + token);
        }
    }

    return linalg::fraction(token);
}

} // namespace

bool parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& A, bool& is_cs)
{
    constexpr size_t kMaxSafeDimension = 63;

    const size_t hash_pos = matrix_str.find('#');
    if (hash_pos == std::string::npos || hash_pos == 0 || hash_pos == matrix_str.length() - 1) {
        std::cerr << "Error: String for the matrix does not include '#' as a separator between dimension and matrix!" << std::endl;
        return false;
    }

    if (matrix_str.find('#', hash_pos + 1) != std::string::npos) {
        std::cerr << "Error: Multiple '#' characters found in matrix string!" << std::endl;
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

    if (n == 0 || n > kMaxSafeDimension) {
        std::cerr << "Error: Safe parser supports dimensions in [1, " << kMaxSafeDimension << "], got " << n << std::endl;
        std::cerr << "Use --unsafe to bypass this validation (not recommended)." << std::endl;
        return false;
    }

    const std::string& values_str = matrix_str.substr(hash_pos + 1);
    std::vector<linalg::fraction> rational_values;
    rational_values.reserve(n * (n + 1) / 2);

    try {
        size_t start = 0;
        size_t comma_pos;

        while (start < values_str.length()) {
            comma_pos = values_str.find(',', start);
            if (comma_pos == std::string::npos) {
                comma_pos = values_str.length();
            }

            rational_values.push_back(parse_fraction_token(values_str, start, comma_pos));

            start = comma_pos + 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: Could not convert matrix values to fraction numbers!" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return false;
    }

    const size_t expected_cs = n / 2;
    const size_t expected_sym = n * (n + 1) / 2;
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

void parse_matrix_string_unsafe(const std::string& matrix_str, linalg::matrix_frc& A, bool& is_cs)
{
    size_t hash_pos = 0;
    while (matrix_str[hash_pos] != '#') {
        ++hash_pos;
    }

    size_t n = 0;
    for (size_t i = 0; i < hash_pos; ++i) {
        n = n * 10 + (matrix_str[i] - '0');
    }

    std::vector<linalg::fraction> rational_values;
    rational_values.reserve(n / 2);

    size_t pos = hash_pos + 1;
    const size_t len = matrix_str.length();

    while (pos < len) {
        int64_t num = 0;
        bool num_negative = false;
        if (matrix_str[pos] == '-') {
            num_negative = true;
            ++pos;
        }
        while (pos < len && matrix_str[pos] != '/' && matrix_str[pos] != ',') {
            num = num * 10 + (matrix_str[pos] - '0');
            ++pos;
        }
        if (num_negative) {
            num = -num;
        }

        if (pos < len && matrix_str[pos] == '/') {
            ++pos;
            int64_t den = 0;
            bool den_negative = false;
            if (matrix_str[pos] == '-') {
                den_negative = true;
                ++pos;
            }
            while (pos < len && matrix_str[pos] != ',') {
                den = den * 10 + (matrix_str[pos] - '0');
                ++pos;
            }
            if (den_negative) {
                num = -num;
            }
            rational_values.push_back(linalg::fraction(static_cast<long long>(num), static_cast<long long>(den)));
        } else {
            rational_values.push_back(linalg::fraction(static_cast<long long>(num), 1LL));
        }
        if (pos < len && matrix_str[pos] == ',') {
            ++pos;
        }
    }

    const size_t expected_cs = n / 2;
    if (rational_values.size() == expected_cs) {
        A = linalg::create_circular_symmetric(n, rational_values);
        is_cs = true;
    } else {
        A = linalg::create_symmetric(n, rational_values);
        is_cs = false;
    }
}

} // namespace matrix_parser
