#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <rational_linalg/matrix_fraction.hpp>
#include <exception>
#include <numeric>

fracessa::fracessa(const rational_linalg::matrix_fraction& matrix, bool is_cs, bool with_candidates, bool exact, bool full_support, bool with_log, int matrix_id)
    : matrix_server_(matrix)
    , dimension_(matrix.rows())
    , is_cs_(is_cs)
    , matrix_id_(matrix_id)
    , conf_with_candidates_(with_candidates)
    , conf_exact_(exact)
    , conf_full_support_(full_support)
    , conf_with_log_(with_log)
    , candidate_()
    , supports_(dimension_, is_cs_)
    , logger_()
{
    if (conf_with_candidates_)
        candidates_.reserve(250 * dimension_);

    if (conf_with_log_) {
        auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            "log/fracessa.log", 20*1024*1024, 5);
        logger_ = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
        logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
        logger_->set_level(spdlog::level::info);
        
        logger_->info("");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");
        logger_->info("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");        
        if (matrix_id >= 0) {
            logger_->info("matrix_id={}", matrix_id);
        }
        
        logger_->info("n={}", dimension_);
        logger_->info("game matrix:\n{}", matrix_server_.get_game_matrix_fraction().to_log_string());
    }
        
    supports_.initialize();

    if (conf_full_support_) {      
        for (const auto& support : supports_.get_supports(dimension_)) {
            search_one_support(support, dimension_);
        }
        if (ess_count_ > 0) 
            return;                
    }

    for (size_t i = 1; i <= (conf_full_support_ ? dimension_-1: dimension_) ; i++) {
        if (conf_with_log_)
            logger_->info("Searching support size {}:", i);

        const bool is_cs_and_coprime = is_cs_ && (std::gcd(i, dimension_) == 1);
        for (const auto& support : supports_.get_supports(i)) {
            search_one_support(support, i, is_cs_and_coprime);
        }
    }
}


void fracessa::search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime)
{
    if (!conf_exact_) 
        if (!find_candidate_double(support, support_size))
            return;

    if (!find_candidate_fraction(support, support_size))
            return;

    candidate_.support_size = support_size;
    candidate_.support = support;
    candidate_.candidate_id++;
    
    if (conf_with_log_)
        logger_->info("Found candidate! Check stability:");

    check_stability();

    if (candidate_.is_ess)
        ess_count_++;

    if (is_cs_and_coprime)
        candidate_.shift_reference = candidate_.candidate_id;
    else
        candidate_.shift_reference = 0;

    if (conf_with_candidates_)
        candidates_.push_back(candidate_);

    if (conf_with_log_) {
        logger_->info("{}", candidate::header());
        logger_->info("{}", candidate_.to_string());
    }
    supports_.remove_supersets(candidate_.support, support_size);

    if (is_cs_and_coprime) {
        for (size_t i = 0; i < dimension_ - 1; i++) {
            candidate_.support = bs64::rot_one_right(candidate_.support, dimension_);
            candidate_.candidate_id++;

            if (candidate_.is_ess)
                ess_count_++;

            if (conf_with_candidates_) {
                fraction first = candidate_.vector(0, 0);
                size_t vec_size = candidate_.vector.rows();
                for (size_t j = 0; j < vec_size - 1; j++) {
                    candidate_.vector(j, 0) = candidate_.vector(j + 1, 0);
                }
                candidate_.vector(vec_size - 1, 0) = first;
                candidate_.extended_support = bs64::rot_one_right(candidate_.extended_support, dimension_);
                candidates_.push_back(candidate_);

                if (conf_with_log_) {
                    logger_->info("{}", candidate::header());
                    logger_->info("{}", candidate_.to_string());
                }
            }
            supports_.remove_supersets(candidate_.support, support_size);
        }
    }
}
