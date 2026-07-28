#pragma once

#include <string>
#include <linalg/matrix_fraction.hpp>

namespace matrix_parser {

/*
 * Parse the shared CLI/Python payload `n#values`. Input requires one
 * separator, a complete decimal dimension in [1, 63], valid exact rationals,
 * and exactly one of these value counts:
 *
 * - floor(n/2) values for the compact circular-symmetric format;
 * - n*(n+1)/2 values for the upper triangle of a symmetric matrix.
 *
 * On success, A contains the expanded square matrix and is_cs identifies the
 * compact format. Callers must ignore both output arguments after `false`.
 */
bool parse_matrix_string(const std::string& matrix_str, linalg::matrix_frc& A, bool& is_cs);

} // namespace matrix_parser
