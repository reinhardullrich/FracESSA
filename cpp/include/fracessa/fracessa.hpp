// Migrated header for modern include path
#include <vector>
#include <string>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <rational_linalg/matrix.hpp>
#include <fracessa/candidate.hpp>
#include <fracessa/bitset64.hpp>
#include <fracessa/supports.hpp>
#include <fracessa/matrix_server.hpp>


class fracessa
{
    public:

        fracessa(const rational_linalg::Matrix<fraction>& matrix, bool is_cs, bool with_candidates = false, bool exact = false, bool full_support = false, bool with_log = false, int matrix_id = -1);

        size_t ess_count_ = 0;
        std::vector<candidate> candidates_;

    private:

        // Matrix server handles all game/augmented/bee matrices and their operations
        MatrixServer matrix_server_;

        size_t dimension_;
        bool is_cs_;
        int matrix_id_;

        bool conf_with_candidates_;
        bool conf_exact_;
        bool conf_full_support_;
        bool conf_with_log_;

        
        candidate candidate_;

        Supports supports_;
        std::vector<bitset64> supports_to_remove_;  // Collect supports for batch removal

        std::shared_ptr<spdlog::logger> logger_;

        void search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime = false);
        
        // Templated function for all types (double, fraction)
        template<typename T>
        bool find_candidate(const bitset64& support, size_t support_size);
        
        // check_stability function for fraction
        void check_stability();

};

