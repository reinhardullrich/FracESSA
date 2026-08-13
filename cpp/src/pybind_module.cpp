#include <pybind11/pybind11.h>

#include <chrono>
#include <cstdint>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <coposit/parsers/matrix_parser.hpp>
#include <fracessa/fracessa.hpp>

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
    std::string support;
    size_t support_size = 0;
    std::string extended_support;
    size_t extended_support_size = 0;
    std::optional<size_t> multiplier;
    bool is_ess = false;
    std::string stability;
    std::string payoff;
    double payoff_dbl = 0.0;
};

struct NativeResult {
    int status = kStatusInternalError;
    std::string error_message;
    size_t candidate_count = 0;
    size_t ess_count = 0;
    std::vector<size_t> candidate_structure;
    std::vector<size_t> ess_structure;
    long long elapsed_ns = 0;
    fracessa::search::safe_fallback safe_fallback = fracessa::search::safe_fallback::none;
    std::vector<NativeCandidate> candidates;
};

std::string strategy_vector_to_string(const std::vector<fracessa::numeric::fraction>& vec)
{
    // Keep exact rational probabilities as text instead of rounding to double.
    std::ostringstream out;
    const size_t n = vec.size();
    for (size_t i = 0; i < n; ++i) {
        out << vec[i].to_string();
        if (i + 1 < n) {
            out << ",";
        }
    }
    return out.str();
}

template<class Analyzer>
void run_analyzer(fracessa::search_method method, coposit::parsers::parsed_matrix matrix, bool include_candidates,
                  bool full_support, bool enable_logging, std::int64_t matrix_id, NativeResult& result)
{
    const auto start = std::chrono::steady_clock::now();
    Analyzer analyzer(method, std::move(matrix), include_candidates, full_support, enable_logging, matrix_id);
    const auto end = std::chrono::steady_clock::now();

    result.elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    result.candidate_count = analyzer.candidate_count_;
    result.ess_count = analyzer.ess_count_;
    result.candidate_structure.assign(analyzer.candidate_structure_.begin(), analyzer.candidate_structure_.end());
    result.ess_structure.assign(analyzer.ess_structure_.begin(), analyzer.ess_structure_.end());
    result.safe_fallback = analyzer.safe_fallback_;

    if (!include_candidates) return;

    result.candidates.reserve(analyzer.candidates_.size());
    for (const auto& candidate : analyzer.candidates_) {
        NativeCandidate row;
        row.candidate_id = candidate.candidate_id;
        row.vector = strategy_vector_to_string(candidate.vector);
        row.support = candidate.support_string();
        row.support_size = candidate.support_size;
        row.extended_support = candidate.extended_support_string();
        row.extended_support_size = candidate.extended_support_size;
        row.multiplier = candidate.multiplier;
        row.is_ess = candidate.is_ess;
        row.stability = candidate.stability;
        row.payoff = candidate.payoff.to_string();
        row.payoff_dbl = candidate.payoff_dbl;
        result.candidates.push_back(std::move(row));
    }
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
        const fracessa::search_method method = fracessa::parse_search_method(method_name);

        coposit::parsers::parsed_matrix parsed_matrix;

        try {
            parsed_matrix = coposit::parsers::matrix_parser::parse(matrix);
        } catch (const std::invalid_argument& e) {
            result.status = kStatusParseError;
            result.error_message = e.what();
            return result;
        }

        std::unique_lock<std::mutex> logging_lock(logging_mutex, std::defer_lock);
        if (enable_logging) {
            logging_lock.lock();
        }

        if (parsed_matrix.matrix.rows() <= fracessa::support::kMaxBitsetDimension)
            run_analyzer<fracessa::analyzer>(method, std::move(parsed_matrix), include_candidates, full_support, enable_logging,
                                             matrix_id, result);
        else
            run_analyzer<fracessa::analyzer_multiword>(method, std::move(parsed_matrix), include_candidates, full_support,
                                                       enable_logging, matrix_id, result);
        result.status = kStatusOk;

        return result;
    } catch (const std::exception& e) {
        result.status = kStatusExecError;
        result.error_message = e.what();
        return result;
    } catch (...) {
        result.status = kStatusInternalError;
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
    out["error_message"] = native.error_message;
    out["candidate_count"] = native.candidate_count;
    out["ess_count"] = native.ess_count;
    py::dict candidate_structure;
    py::dict ess_structure;
    for (size_t support_size = 1; support_size < native.candidate_structure.size(); ++support_size) {
        if (native.candidate_structure[support_size])
            candidate_structure[py::int_(support_size)] = native.candidate_structure[support_size];
        if (native.ess_structure[support_size])
            ess_structure[py::int_(support_size)] = native.ess_structure[support_size];
    }
    out["candidate_structure"] = std::move(candidate_structure);
    out["ess_structure"] = std::move(ess_structure);
    out["elapsed_ns"] = native.elapsed_ns;
    const auto fallback = fracessa::search::safe_fallback_name(native.safe_fallback);
    if (fallback.empty())
        out["safe_fallback"] = py::none();
    else
        out["safe_fallback"] = std::string(fallback);

    py::list candidates;
    for (const auto& c : native.candidates) {
        py::dict row;
        row["candidate_id"] = c.candidate_id;
        row["vector"] = c.vector;
        row["support"] = py::int_(py::str(c.support));
        row["support_size"] = c.support_size;
        row["extended_support"] = py::int_(py::str(c.extended_support));
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

fracessa::search::safe_fallback classify_safe_fallback_impl(const std::string& matrix)
{
    const auto parsed_matrix = coposit::parsers::matrix_parser::parse(matrix);
    fracessa::search::fast_candidate_filter fast_filter(parsed_matrix.matrix.rows());
    fast_filter.convert_game_matrix(parsed_matrix.matrix);
    return fast_filter.safe_fallback_reason();
}

py::object classify_safe_fallback(const std::string& matrix)
{
    fracessa::search::safe_fallback fallback;
    {
        py::gil_scoped_release release;
        fallback = classify_safe_fallback_impl(matrix);
    }
    const auto name = fracessa::search::safe_fallback_name(fallback);
    if (name.empty()) return py::none();
    return py::str(name);
}

} // namespace

PYBIND11_MODULE(fracessa_core, m)
{
    m.doc() = "Native pybind11 interface for FracESSA";

    m.attr("STATUS_OK") = kStatusOk;
    m.attr("STATUS_PARSE_ERROR") = kStatusParseError;
    m.attr("STATUS_EXEC_ERROR") = kStatusExecError;
    m.attr("STATUS_INTERNAL_ERROR") = kStatusInternalError;

    m.def("classify_safe_fallback", &classify_safe_fallback, py::arg("matrix"),
          "Return the whole-matrix fast-to-safe fallback reason, or None.");

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
- status, error_message
- candidate_count, ess_count, candidate_structure, ess_structure
- elapsed_ns, safe_fallback
- candidates (list of dict rows)
)doc"
    );
}
