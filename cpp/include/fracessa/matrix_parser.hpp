#pragma once

#include <string>
#include <linalg/matrix_fraction.hpp>

namespace matrix_parser {

// Safe parser for CLI matrix payload `n#values`.
// Accepts:
// - circular-symmetric compact values: floor(n/2)
// - upper-triangular symmetric values: n*(n+1)/2
// Returns false on malformed input.
bool parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& A, bool& is_cs);

// Unsafe parser variant: minimal validation and branch-light scanning.
// Intended for trusted high-throughput input paths.
void parse_matrix_string_unsafe(const std::string& matrix_str, linalg::matrix_frc& A, bool& is_cs);

} // namespace matrix_parser

