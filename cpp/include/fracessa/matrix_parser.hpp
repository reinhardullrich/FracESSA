#pragma once

#include <string>
#include <linalg/matrix_fraction.hpp>

namespace matrix_parser {

// Safe parser for CLI matrix payload `n#values`.
// Accepts:
// - circular-symmetric compact values: floor(n/2)
// - upper-triangular symmetric values: n*(n+1)/2
// Throws std::invalid_argument on malformed input.
void parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& A, bool& is_cs);

} // namespace matrix_parser
