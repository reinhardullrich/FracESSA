#include <pybind11/pybind11.h>

#include <chrono>
#include <cstdint>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fracessa/fracessa.hpp>
#include <fracessa/matrix_parser.hpp>

namespace py = pybind11;

namespace {

// These integer values are part of the Python wrapper's public status contract.
constexpr int kStatusOk = 0;
constexpr int kStatusParseError = 1;
constexpr int kStatusExecError = 4;
constexpr int kStatusInternalError = 255;

// ponytail: one process-wide lock is enough; use per-run log files if logging
// throughput matters.
static std::mutex logging_mutex;

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
    size_t ess_count = 0;
    long long elapsed_ns = 0;
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

NativeResult compute_matrix_impl(
    const std::string& method_name,
    const std::string& matrix,
    bool include_candidates,
    bool full_support,
    bool enable_logging,
    std::int64_t matrix_id)
{
    NativeResult result;

    try {
        const search_method method = parse_search_method(method_name);

        linalg::matrix_frc parsed_matrix;
        bool is_cs = false;

        try {
            matrix_parser::parse_matrix_string(matrix, parsed_matrix, is_cs);
        } catch (const std::invalid_argument& e) {
            result.status = kStatusParseError;
            result.success = false;
            result.error_message = e.what();
            return result;
        }

        std::unique_lock<std::mutex> logging_lock(logging_mutex, std::defer_lock);
        if (enable_logging) {
            logging_lock.lock();
        }

        const auto start = std::chrono::steady_clock::now();
        ::fracessa analyzer(method, parsed_matrix, is_cs, include_candidates, full_support, enable_logging, matrix_id);
        const auto end = std::chrono::steady_clock::now();

        result.elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        result.ess_count = analyzer.ess_count_;
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
    const std::string& method,
    const std::string& matrix,
    bool include_candidates,
    bool full_support,
    bool enable_logging,
    std::int64_t matrix_id)
{
    NativeResult native;
    {
        // compute_matrix_impl touches no Python objects and may run for hours.
        py::gil_scoped_release release;
        native = compute_matrix_impl(method, matrix, include_candidates, full_support, enable_logging, matrix_id);
    }

    // The release object above has restored the GIL; Python allocation is safe.
    py::dict out;
    out["status"] = native.status;
    out["success"] = native.success;
    out["error_message"] = native.error_message;
    out["ess_count"] = native.ess_count;
    out["elapsed_ns"] = native.elapsed_ns;

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
    m.attr("STATUS_EXEC_ERROR") = kStatusExecError;
    m.attr("STATUS_INTERNAL_ERROR") = kStatusInternalError;

    m.def(
        "compute_matrix",
        &compute_matrix,
        py::arg("method"),
        py::arg("matrix"),
        py::arg("include_candidates") = false,
        py::arg("full_support") = false,
        py::arg("enable_logging") = false,
        py::arg("matrix_id") = std::int64_t{-1},
        R"doc(
Compute one matrix with native C++ core and return structured results.

method: fast or safe; required before matrix.

Returns a dict with keys:
- status, success, error_message
- ess_count, elapsed_ns
- candidates (list of dict rows)
)doc"
    );
}
