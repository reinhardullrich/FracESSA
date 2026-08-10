#pragma once

#include <string>

#include <linalg/matrix_fraction.hpp>

namespace matrix_parser {

/// Parse CLI/Python payload `n#values` into an exact symmetric matrix. The value list may contain either floor(n/2) compact
/// circular-symmetric entries or n*(n+1)/2 upper-triangular entries. Malformed input throws std::invalid_argument.
void parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& matrix, bool& is_circular_symmetric);

} // namespace matrix_parser
