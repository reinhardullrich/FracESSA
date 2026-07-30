#include <iostream>
#include <string>
#include <chrono>
#include <cstdint>
#include <stdexcept>

#include <fracessa/fracessa.hpp>
#include <fracessa/matrix_parser.hpp>
#include <argparse/argparse.hpp>

int main(int argc, char *argv[])
{
    /*
     * This file is only the command-line adapter: parse options and one matrix,
     * run the native analyzer, then print its stable line-oriented result.
     */
    argparse::ArgumentParser program("fracessa", "3.0.0");
    program.add_description("FRACESSA - Fractional ESS Analyzer");

    program.add_argument("-c", "--candidates").help("include candidates").implicit_value(true).default_value(false);
    program.add_argument("-l", "--log").help("output log file").implicit_value(true).default_value(false);
    program.add_argument("-e", "--exact").help("exact only").implicit_value(true).default_value(false);
    program.add_argument("-f", "--fullsupport").help("search full support directly").implicit_value(true).default_value(false);
    program.add_argument("-t", "--timing").help("output computation time in nanoseconds").implicit_value(true).default_value(false);
    program.add_argument("-m", "--matrixid").help("optional matrix ID").scan<'i', std::int64_t>().default_value(std::int64_t{-1});
    program.add_argument("-u", "--unsafe").help("use heuristic candidate search").implicit_value(true).default_value(false);
    program.add_argument("matrix").help("the matrix to compute");

    try { program.parse_args(argc, argv); }
    catch (const std::exception& err) { std::cerr << err.what() << std::endl << program; return EXIT_FAILURE; }

    const auto& matrix_str = program.get<std::string>("matrix");
    const auto candidates = program.get<bool>("--candidates");
    const auto logger = program.get<bool>("--log");
    const auto exact = program.get<bool>("--exact");
    const auto fullsupport = program.get<bool>("--fullsupport");
    const auto timing = program.get<bool>("--timing");
    const auto matrix_id = program.get<std::int64_t>("--matrixid");
    const auto unsafe = program.get<bool>("--unsafe");

    if (exact && unsafe) {
        std::cerr << "Error: --exact and --unsafe cannot be used together." << std::endl;
        return EXIT_FAILURE;
    }

    bool is_cs;
    linalg::matrix_frc A;
    try { matrix_parser::parse_matrix_string(matrix_str, A, is_cs); }
    catch (const std::invalid_argument& err) {
        std::cerr << "Error: " << err.what() << std::endl;
        return EXIT_FAILURE;
    }

    if (unsafe) {
        std::cerr
            << "Warning: --unsafe uses heuristic candidate search and can miss "
            << "exact candidates and ESS results." << std::endl;
    }
    try {
        auto start_time = std::chrono::steady_clock::now();
        ::fracessa x(A, is_cs, candidates, exact, fullsupport, logger, matrix_id, unsafe);
        auto end_time = std::chrono::steady_clock::now();

        // Consumers expect ESS count first, optional timing second, then candidate CSV.
        std::cout << x.ess_count_ << std::endl;
        if (timing) {
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
            std::cout << duration.count() << std::endl;
        }

        if (candidates) {
            std::cout << candidate::header() << std::endl;
            for (auto& c : x.candidates_) std::cout << c.to_string() << std::endl;
        }
    } catch (const std::exception& err) {
        std::cerr << "Error: " << err.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
