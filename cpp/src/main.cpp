#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>

#include <argparse/argparse.hpp>
#include <coposit/parsers/matrix_parser.hpp>
#include <fracessa/fracessa.hpp>

namespace {

constexpr int kStatusOk = 0;
constexpr int kStatusParseError = 1;
constexpr int kStatusExecError = 4;
constexpr int kStatusInternalError = 255;

void write_json_string(std::ostream& output, std::string_view value)
{
    static constexpr char hex[] = "0123456789abcdef";
    output << '"';
    for (const unsigned char character : value) {
        switch (character) {
            case '"': output << "\\\""; break;
            case '\\': output << "\\\\"; break;
            case '\b': output << "\\b"; break;
            case '\f': output << "\\f"; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default:
                if (character < 0x20)
                    output << "\\u00" << hex[character >> 4] << hex[character & 0x0f];
                else
                    output << static_cast<char>(character);
        }
    }
    output << '"';
}

template<class Structure>
void write_structure(std::ostream& output, const Structure& structure)
{
    output << '{';
    bool first = true;
    for (size_t support_size = 1; support_size < structure.size(); ++support_size) {
        if (!structure[support_size]) continue;
        if (!first) output << ',';
        output << '"' << support_size << "\":" << structure[support_size];
        first = false;
    }
    output << '}';
}

template<class Structure>
void write_summary(
    std::int64_t matrix_id,
    int status,
    size_t candidate_count,
    size_t ess_count,
    const Structure& candidate_structure,
    const Structure& ess_structure,
    long long elapsed_ns,
    fracessa::search::safe_fallback safe_fallback,
    std::string_view error_message)
{
    std::cout << R"({"run_id":null,"matrix_id":)" << matrix_id
              << R"(,"status":)" << status
              << R"(,"candidate_count":)" << candidate_count
              << R"(,"ess_count":)" << ess_count
              << R"(,"candidate_structure":)";
    write_structure(std::cout, candidate_structure);
    std::cout << R"(,"ess_structure":)";
    write_structure(std::cout, ess_structure);
    std::cout << R"(,"elapsed_ns":)" << elapsed_ns << R"(,"safe_fallback":)";
    const auto fallback_name = fracessa::search::safe_fallback_name(safe_fallback);
    if (fallback_name.empty())
        std::cout << "null";
    else
        write_json_string(std::cout, fallback_name);
    std::cout << R"(,"error_message":)";
    write_json_string(std::cout, error_message);
    std::cout << '}' << std::endl;
}

void write_error(std::int64_t matrix_id, int status, std::string_view message)
{
    const std::array<size_t, fracessa::support::kMaxBitsetDimension + 1> empty_structure{};
    write_summary(
        matrix_id,
        status,
        0,
        0,
        empty_structure,
        empty_structure,
        0,
        fracessa::search::safe_fallback::none,
        message);
}

template<class Analyzer>
void run_analysis(fracessa::search_method method, coposit::parsers::parsed_matrix matrix, bool include_candidates,
                  bool full_support, bool enable_logging, std::int64_t matrix_id)
{
    const auto start_time = std::chrono::steady_clock::now();
    Analyzer analyzer(method, std::move(matrix), include_candidates, full_support, enable_logging, matrix_id);
    const auto end_time = std::chrono::steady_clock::now();
    const auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

    write_summary(matrix_id, kStatusOk, analyzer.candidate_count_, analyzer.ess_count_, analyzer.candidate_structure_,
                  analyzer.ess_structure_, elapsed_ns, analyzer.safe_fallback_, "");

    if (include_candidates) {
        using candidate_type = typename Analyzer::candidate_type;
        std::cout << candidate_type::header() << std::endl;
        for (const auto& row : analyzer.candidates_) std::cout << row.to_string() << std::endl;
    }
}

bool is_inline_matrix(std::string_view argument) noexcept
{
    const size_t separator = argument.find('#');
    if (separator == 0 || separator == std::string_view::npos) return false;
    for (size_t position = 0; position < separator; ++position) {
        if (argument[position] < '0' || argument[position] > '9') return false;
    }
    return true;
}

} // namespace

int main(int argc, char *argv[])
{
    /*
     * This file is only the command-line adapter: parse one matrix, run the
     * native analyzer, then print one JSON summary and optional candidate CSV.
     */
    argparse::ArgumentParser program("fracessa", FRACESSA_VERSION);
    program.add_description("FRACESSA - Fractional ESS Analyzer");

    program.add_argument("-c", "--candidates").help("include candidate CSV").implicit_value(true).default_value(false);
    program.add_argument("-l", "--log").help("output log file").implicit_value(true).default_value(false);
    program.add_argument("-f", "--fullsupport").help("search full support directly").implicit_value(true).default_value(false);
    program.add_argument("-m", "--matrixid").help("optional matrix ID").scan<'i', std::int64_t>().default_value(std::int64_t{-1});
    program.add_argument("method").help("candidate search method: fast or safe");
    program.add_argument("matrix").help("inline dimension#values matrix or matrix file");

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& error) {
        write_error(-1, kStatusParseError, error.what());
        std::cerr << error.what() << std::endl << program;
        return EXIT_FAILURE;
    }

    const auto& method_name = program.get<std::string>("method");
    const auto& matrix_argument = program.get<std::string>("matrix");
    const bool include_candidates = program.get<bool>("--candidates");
    const bool enable_logging = program.get<bool>("--log");
    const bool full_support = program.get<bool>("--fullsupport");
    const std::int64_t matrix_id = program.get<std::int64_t>("--matrixid");

    fracessa::search_method method;
    try {
        method = fracessa::parse_search_method(method_name);
    } catch (const std::invalid_argument& error) {
        write_error(matrix_id, kStatusExecError, error.what());
        std::cerr << "Error: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }

    coposit::parsers::parsed_matrix matrix;
    try {
        matrix = is_inline_matrix(matrix_argument) ? coposit::parsers::matrix_parser::parse(matrix_argument)
                                                   : coposit::parsers::matrix_parser::parse_file(matrix_argument);
    } catch (const std::invalid_argument& error) {
        write_error(matrix_id, kStatusParseError, error.what());
        std::cerr << "Error: " << error.what() << std::endl;
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        write_error(matrix_id, kStatusExecError, error.what());
        std::cerr << "Error: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }

    try {
        if (matrix.matrix.rows() <= fracessa::support::kMaxBitsetDimension)
            run_analysis<fracessa::analyzer>(method, std::move(matrix), include_candidates, full_support, enable_logging, matrix_id);
        else
            run_analysis<fracessa::analyzer_multiword>(method, std::move(matrix), include_candidates, full_support, enable_logging,
                                                       matrix_id);
    } catch (const std::exception& error) {
        write_error(matrix_id, kStatusExecError, error.what());
        std::cerr << "Error: " << error.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        constexpr std::string_view message = "Unknown internal error";
        write_error(matrix_id, kStatusInternalError, message);
        std::cerr << "Error: " << message << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
