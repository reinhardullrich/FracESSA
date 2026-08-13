#pragma once

#include <cstddef>
#include <sstream>
#include <string>

#include <coposit/matrix_integer.hpp>
#include <linalg/integer.hpp>

namespace linalg {

using matrix_int = coposit::matrix_integer;

inline std::string to_pretty_string(const matrix_int& matrix)
{
    std::stringstream stream;
    for (size_t row = 0; row < matrix.rows(); ++row) {
        stream << "  " << row << ": [";
        for (size_t column = 0; column < matrix.cols(); ++column) {
            if (column != 0) stream << ", ";
            stream << coposit::integer(matrix(row, column)).to_string();
        }
        stream << ']';
        if (row + 1 != matrix.rows()) stream << '\n';
    }
    return stream.str();
}

} // namespace linalg
