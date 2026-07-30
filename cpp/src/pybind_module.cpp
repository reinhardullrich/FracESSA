#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <chrono>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <fracessa/fracessa.hpp>
#include <fracessa/matrix_parser.hpp>

namespace py = pybind11;

namespace {

// These integer values are part of the Python wrapper's public status contract.
constexpr int kStatusOk = 0;
constexpr int kStatusParseError = 1;
constexpr int kStatusDimensionOutOfRange = 2;
constexpr int kStatusInvalidValueCount = 3;
constexpr int kStatusExecError = 4;
constexpr int kStatusInternalError = 255;

/*
 * Python-free transport objects let parsing, analysis, and candidate copying
 * run while the GIL is released. Python dictionaries are created only after
 * native work has finished.
 */
struct NativeCandidate {
    size_t candidate_id = 0;
    std::string vector;
    bitset64 support = 0;
    size_t support_size = 0;
    bitset64 extended_support = 0;
    size_t extended_support_size = 0;
    std::optional<size_t> multiplier;
    bool is_ess = false;
    std::string stability;
    std::string payoff;
    double payoff_dbl = 0.0;
};

struct NativeResult {
    int status = kStatusInternalError;
    bool success = false;
    std::string error_message;
    int ess_count = 0;
    long long elapsed_us = 0;
    std::vector<NativeCandidate> candidates;
};

std::string strategy_vector_to_string(const linalg::matrix_frc& vec)
{
    // Keep exact rational probabilities as text instead of rounding to double.
    std::ostringstream out;
    const size_t n = vec.rows();
    for (size_t i = 0; i < n; ++i) {
        out << vec(i, 0).to_string();
        if (i + 1 < n) {
            out << ",";
        }
    }
    return out.str();
}

int infer_parse_status(const std::string& matrix_str)
{
    /*
     * The shared parser currently returns only bool. This best-effort
     * second pass recognizes some structural and dimension failures; every
     * other parse failure is reported as an invalid value count.
     */
    constexpr size_t kMaxSafeDimension = 63;

    const size_t hash_pos = matrix_str.find('#');
    if (hash_pos == std::string::npos || hash_pos == 0 || hash_pos == matrix_str.length() - 1) {
        return kStatusParseError;
    }
    if (matrix_str.find('#', hash_pos + 1) != std::string::npos) {
        return kStatusParseError;
    }

    try {
        const size_t n = std::stoull(matrix_str.substr(0, hash_pos));
        if (n == 0 || n > kMaxSafeDimension) {
            return kStatusDimensionOutOfRange;
        }
    } catch (const std::exception&) {
        return kStatusParseError;
    }

    return kStatusInvalidValueCount;
}

NativeResult compute_matrix_impl(
    const std::string& matrix,
    bool include_candidates,
    bool exact,
    bool full_support,
    bool enable_logging,
    int matrix_id)
{
    NativeResult result;

    try {
        linalg::matrix_frc parsed_matrix;
        bool is_cs = false;

        if (!matrix_parser::parse_matrix_string(matrix, parsed_matrix, is_cs)) {
            result.status = infer_parse_status(matrix);
            result.success = false;
            result.error_message = "Failed to parse matrix input";
            return result;
        }

        // Match CLI timing: measure only analyzer construction and search.
        const auto start = std::chrono::high_resolution_clock::now();
        ::fracessa analyzer(parsed_matrix, is_cs, include_candidates, exact, full_support, enable_logging, matrix_id);
        const auto end = std::chrono::high_resolution_clock::now();

        result.elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        result.ess_count = static_cast<int>(analyzer.ess_count_);
        result.status = kStatusOk;
        result.success = true;

        if (include_candidates) {
            // Copy all C++ data before the GIL is reacquired below.
            result.candidates.reserve(analyzer.candidates_.size());
            for (const auto& c : analyzer.candidates_) {
                NativeCandidate row;
                row.candidate_id = c.candidate_id;
                row.vector = strategy_vector_to_string(c.vector);
                row.support = c.support;
                row.support_size = c.support_size;
                row.extended_support = c.extended_support;
                row.extended_support_size = c.extended_support_size;
                row.multiplier = c.multiplier;
                row.is_ess = c.is_ess;
                row.stability = c.stability;
                row.payoff = c.payoff.to_string();
                row.payoff_dbl = c.payoff_dbl;
                result.candidates.push_back(std::move(row));
            }
        }

        return result;
    } catch (const std::exception& e) {
        result.status = kStatusExecError;
        result.success = false;
        result.error_message = e.what();
        return result;
    } catch (...) {
        result.status = kStatusInternalError;
        result.success = false;
        result.error_message = "Unknown internal error";
        return result;
    }
}

py::dict compute_matrix(
    const std::string& matrix,
    bool include_candidates,
    bool exact,
    bool full_support,
    bool enable_logging,
    int matrix_id)
{
    NativeResult native;
    {
        // compute_matrix_impl touches no Python objects and may run for hours.
        py::gil_scoped_release release;
        native = compute_matrix_impl(
            matrix, include_candidates, exact, full_support,
            enable_logging, matrix_id);
    }

    // The release object above has restored the GIL; Python allocation is safe.
    py::dict out;
    out["status"] = native.status;
    out["success"] = native.success;
    out["error_message"] = native.error_message;
    out["ess_count"] = native.ess_count;
    out["elapsed_us"] = native.elapsed_us;

    py::list candidates;
    for (const auto& c : native.candidates) {
        py::dict row;
        row["candidate_id"] = c.candidate_id;
        row["vector"] = c.vector;
        row["support"] = c.support;
        row["support_size"] = c.support_size;
        row["extended_support"] = c.extended_support;
        row["extended_support_size"] = c.extended_support_size;
        if (c.multiplier)
            row["multiplier"] = *c.multiplier;
        else
            row["multiplier"] = py::none();
        row["is_ess"] = c.is_ess;
        row["stability"] = c.stability;
        row["payoff"] = c.payoff;
        row["payoff_dbl"] = c.payoff_dbl;
        candidates.append(std::move(row));
    }
    out["candidates"] = std::move(candidates);
    return out;
}

} // namespace

PYBIND11_MODULE(fracessa_core, m)
{
    m.doc() = "Native pybind11 interface for FracESSA";

    m.attr("STATUS_OK") = kStatusOk;
    m.attr("STATUS_PARSE_ERROR") = kStatusParseError;
    m.attr("STATUS_DIMENSION_OUT_OF_RANGE") = kStatusDimensionOutOfRange;
    m.attr("STATUS_INVALID_VALUE_COUNT") = kStatusInvalidValueCount;
    m.attr("STATUS_EXEC_ERROR") = kStatusExecError;
    m.attr("STATUS_INTERNAL_ERROR") = kStatusInternalError;

    m.def(
        "compute_matrix",
        &compute_matrix,
        py::arg("matrix"),
        py::arg("include_candidates") = false,
        py::arg("exact") = false,
        py::arg("full_support") = false,
        py::arg("enable_logging") = false,
        py::arg("matrix_id") = -1,
        R"doc(
Compute one matrix with native C++ core and return structured results.

Returns a dict with keys:
- status, success, error_message
- ess_count, elapsed_us
- candidates (list of dict rows)
)doc"
    );
}
