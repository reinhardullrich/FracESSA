#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

#include <fracessa/fracessa.hpp>
#include <fracessa/matrix_parser.hpp>
#include <argparse/argparse.hpp>

int main(int argc, char *argv[])
{
    // CLI surface mirrors core search toggles used in batch verification and CI.
    argparse::ArgumentParser program("fracessa", "3.0.0");
    program.add_description("FRACESSA - Fractional ESS Analyzer");

    program.add_argument("-c", "--candidates").help("include candidates").implicit_value(true).default_value(false);
    program.add_argument("-l", "--log").help("output log file").implicit_value(true).default_value(false);
    program.add_argument("-e", "--exact").help("exact only").implicit_value(true).default_value(false);
    program.add_argument("-f", "--fullsupport").help("search full support directly").implicit_value(true).default_value(false);
    program.add_argument("-t", "--timing").help("output computation time").implicit_value(true).default_value(false);
    program.add_argument("-m", "--matrixid").help("optional matrix ID").scan<'i', int>().default_value(-1);
    program.add_argument("-u", "--unsafe").help("unsafe parsing").implicit_value(true).default_value(false);
    program.add_argument("matrix").help("the matrix to compute");

    try { program.parse_args(argc, argv); }
    catch (const std::exception& err) { std::cerr << err.what() << std::endl << program; return EXIT_FAILURE; }

    const auto& matrix_str = program.get<std::string>("matrix");
    const auto candidates = program.get<bool>("--candidates");
    const auto logger = program.get<bool>("--log");
    const auto exact = program.get<bool>("--exact");
    const auto fullsupport = program.get<bool>("--fullsupport");
    const auto timing = program.get<bool>("--timing");
    const auto matrix_id = program.get<int>("--matrixid");
    const auto unsafe = program.get<bool>("--unsafe");

    bool is_cs;
    linalg::matrix_frc A;
    if (unsafe) matrix_parser::parse_matrix_string_unsafe(matrix_str, A, is_cs);
    else if (!matrix_parser::parse_matrix_string(matrix_str, A, is_cs)) return EXIT_FAILURE;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    ::fracessa x(A, is_cs, candidates, exact, fullsupport, logger, matrix_id);
    auto end_time = std::chrono::high_resolution_clock::now();
    
    std::cout << x.ess_count_ << std::endl;
    if (timing) {
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        std::cout << std::fixed << std::setprecision(6) << duration.count() / 1000000.0 << std::endl;
    }

    if (candidates) {
        std::cout << candidate::header() << std::endl;
        for (auto& c : x.candidates_) std::cout << c.to_string() << std::endl;
    }

    return EXIT_SUCCESS;
}
