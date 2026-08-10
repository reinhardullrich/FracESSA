#pragma once

#include <string>

#include <linalg/matrix_fraction.hpp>

namespace matrix_parser {

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
void parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& matrix, bool& is_circular_symmetric);

} // namespace matrix_parser
