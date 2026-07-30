#include <iostream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <string>
#include <chrono>
#include <iomanip>

#include "../include/EssFinder.hpp"
#include "../include/Matrix.hpp"
#include "../include/argtable3.h"

struct arg_lit *candidates, *help, *version, *logger, *exact, *fullsupport, *timing;
struct arg_str *matrixinput;
struct arg_end *end;

int main(int argc, char *argv[])
{

    void *argtable[] = {
        help            = arg_litn("h", "help", 0, 1, "display this help and exit"),
        //version         = arg_litn(NULL, "version", 0, 1, "display version info and exit"),
        candidates         = arg_litn("c", "candidates", 0, 1, "include the found candidates for ESS/solutions in the output"),
        logger          = arg_litn("l", "log", 0, 1, "output a detailed log file named 'essfinder_log.txt' in the directory of the program, for diagnostic of learning purposes only"),
        exact           = arg_litn("e", "exact", 0, 1, "only uses rational numbers, for matrices with extreme differences in the input, is much much slower!"),
        fullsupport     = arg_litn("f", "fullsupport", 0, 1, "searches the full support directly after searching support size one. Enable if you expect the matrix to have exactly one ess in the interior of the simplex!"),
        timing          = arg_litn("t", "timing", 0, 1, "output execution time in seconds on the second line"),
        matrixinput     = arg_strn(NULL, NULL, "<matrix>", 1, 1, "the matrix to compute"),
        end     = arg_end(50),
    };

    char progname[] = "REF.exe";

    int nerrors;
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0) {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        printf("Demonstrate command-line parsing in argtable3.\n\n");
        arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        exit(EXIT_SUCCESS);
    }
    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout, end, progname);
        printf("Try '%s --help' for more information.\n", progname);
        exit(EXIT_FAILURE);
    }

    Matrix<rational> A = Matrix<rational>::create_from_cli_string(matrixinput->sval[0]);
    
    // Start timing
    const auto start_time = std::chrono::high_resolution_clock::now();
    
    EssFinder x = EssFinder(A, candidates->count, exact->count, fullsupport->count, logger->count);
    
    // End timing
    const auto end_time = std::chrono::high_resolution_clock::now();
    
    // Calculate elapsed time in seconds
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    const double elapsed_seconds = duration.count() / 1000000.0;

    std::cout << x.ess_count << std::endl;
    
    // Output timing on second line if -t flag is present
    if (timing->count > 0) {
        std::cout << std::fixed << std::setprecision(6) << elapsed_seconds << std::endl;
    }

    if (candidates->count) {
        std::cout << Candidate::header() << std::endl;
        for (auto c: x.candidates) {
            std::cout << c.to_string() << std::endl;
        }
    }

}
