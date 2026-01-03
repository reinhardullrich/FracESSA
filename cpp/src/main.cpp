#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <cstdint>

#include <fracessa/fracessa.hpp>
#include <rational_linalg/matrix_fraction.hpp>
#include <argparse/argparse.hpp>

// Helper function to parse matrix string format: "n#values"
bool parse_matrix_string(const std::string& matrix_str, rational_linalg::matrix_fraction& A, bool& is_cs)
{
    const size_t hash_pos = matrix_str.find('#');
    if (hash_pos == std::string::npos || hash_pos == 0 || hash_pos == matrix_str.length() - 1) {
        std::cerr << "Error: String for the matrix does not include '#' as a separator between dimension and matrix!" << std::endl;
        return false;
    }
    
    if (matrix_str.find('#', hash_pos + 1) != std::string::npos) {
        std::cerr << "Error: Multiple '#' characters found in matrix string!" << std::endl;
        return false;
    }
    
    size_t n;
    try {
        n = std::stoull(matrix_str.substr(0, hash_pos));
    } catch (const std::exception& e) {
        std::cerr << "Error: The given dimension could not be converted into an integer number!" << std::endl;
        return false;
    }
    
    const std::string& values_str = matrix_str.substr(hash_pos + 1);
    std::vector<fraction> rational_values;
    rational_values.reserve(n / 2);
    
    try {
        size_t start = 0;
        size_t comma_pos;
        
        while (start < values_str.length()) {
            comma_pos = values_str.find(',', start);
            if (comma_pos == std::string::npos) {
                comma_pos = values_str.length();
            }
            
            const size_t slash_pos = values_str.find('/', start);
            if (slash_pos != std::string::npos && slash_pos < comma_pos) {
                int64_t num = std::stoll(values_str.substr(start, slash_pos - start));
                int64_t den = std::stoll(values_str.substr(slash_pos + 1, comma_pos - slash_pos - 1));
                rational_values.push_back(fraction(static_cast<long>(num), static_cast<long>(den)));
            } else {
                int64_t num = std::stoll(values_str.substr(start, comma_pos - start));
                rational_values.push_back(fraction(static_cast<long>(num), 1L));
            }
            
            start = comma_pos + 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: Could not convert matrix values to fraction numbers!" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return false;
    }
    
    const size_t expected_cs = n / 2;
    const size_t expected_sym = n * (n + 1) / 2;
    const size_t actual_size = rational_values.size();
    
    if (actual_size == expected_cs) {
        A = rational_linalg::create_circular_symmetric(n, rational_values);
        is_cs = true;
    } else if (actual_size == expected_sym) {
        A = rational_linalg::create_symmetric(n, rational_values);
        is_cs = false;
    } else {
        std::cerr << "Error: Expected " << expected_cs << " (CS) or " << expected_sym << " (Sym) values, got " << actual_size << std::endl;
        return false;
    }
    
    return true;
}

void parse_matrix_string_unsafe(const std::string& matrix_str, rational_linalg::matrix_fraction& A, bool& is_cs)
{
    size_t hash_pos = 0;
    while (matrix_str[hash_pos] != '#') ++hash_pos;
    
    size_t n = 0;
    for (size_t i = 0; i < hash_pos; ++i) n = n * 10 + (matrix_str[i] - '0');
    
    std::vector<fraction> rational_values;
    rational_values.reserve(n / 2);
    
    size_t pos = hash_pos + 1;
    const size_t len = matrix_str.length();
    
    while (pos < len) {
        int64_t num = 0;
        bool num_negative = false;
        if (matrix_str[pos] == '-') { num_negative = true; ++pos; }
        while (pos < len && matrix_str[pos] != '/' && matrix_str[pos] != ',') {
            num = num * 10 + (matrix_str[pos] - '0');
            ++pos;
        }
        if (num_negative) num = -num;
        
        if (pos < len && matrix_str[pos] == '/') {
            ++pos;
            int64_t den = 0;
            bool den_negative = false;
            if (matrix_str[pos] == '-') { den_negative = true; ++pos; }
            while (pos < len && matrix_str[pos] != ',') {
                den = den * 10 + (matrix_str[pos] - '0');
                ++pos;
            }
            if (den_negative) den = -den;
            rational_values.push_back(fraction(static_cast<long>(num), static_cast<long>(den)));
        } else {
            rational_values.push_back(fraction(static_cast<long>(num), 1L));
        }
        if (pos < len && matrix_str[pos] == ',') ++pos;
    }
    
    const size_t expected_cs = n / 2;
    if (rational_values.size() == expected_cs) {
        A = rational_linalg::create_circular_symmetric(n, rational_values);
        is_cs = true;
    } else {
        A = rational_linalg::create_symmetric(n, rational_values);
        is_cs = false;
    }
}

int main(int argc, char *argv[])
{
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
    rational_linalg::matrix_fraction A;
    if (unsafe) parse_matrix_string_unsafe(matrix_str, A, is_cs);
    else if (!parse_matrix_string(matrix_str, A, is_cs)) return EXIT_FAILURE;
    
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
