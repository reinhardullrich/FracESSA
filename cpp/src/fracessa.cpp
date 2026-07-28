#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/matrix_fraction.hpp>
#include <numeric>

/*
 * End-to-end search for one game.
 *
 * A support first has to produce an exact symmetric Nash equilibrium. Such an
 * equilibrium is called a candidate here. Only then is its stronger ESS
 * stability condition checked. Supports are visited by increasing size so that
 * every exact equilibrium can prune impossible strict supersets before those
 * larger supports are reached.
 */

fracessa::fracessa(const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates, bool exact, bool full_support, bool with_log, int matrix_id)
    : matrix_server_(matrix)
    , dimension_(matrix.rows())
    , is_cs_(is_cs)
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
        logger_->info("game matrix:\n{}", matrix_server_.get_game_matrix_frc().to_log_string());
    }

    // Normalization is done once, not once per support. If binary64 cannot
    // represent the normalized game safely enough for this heuristic, use the
    // exact candidate path for the complete search.
    if (!conf_exact_ && !matrix_server_.initialize_unsafe_filter())
        conf_exact_ = true;
        
    // --fullsupport is a cheap first attempt: test its one mask before spending
    // time and memory constructing every support. If it is not an ESS, continue
    // with all smaller supports as the documented fallback.
    if (conf_full_support_) {
        const bitset64 full_support_mask = bs64::set_all_n_bits(dimension_);
        search_one_support(full_support_mask, dimension_);
        if (ess_count_ > 0)
            return;
    }

    // Build all support buckets only when the normal or fallback search needs them.
    supports_.initialize();

    // Increasing cardinality is essential for valid superset pruning.
    const size_t max_support_size = conf_full_support_ ? dimension_ - 1 : dimension_;
    for (size_t i = 1; i <= max_support_size; i++) {
        if (conf_with_log_)
            logger_->info("Searching support size {}:", i);

        // In this case every rotational orbit has n distinct members, so only
        // its canonical support was stored and the other members can be copied.
        const bool is_cs_and_coprime = is_cs_ && (std::gcd(i, dimension_) == 1);
        for (const auto& support : supports_.get_supports(i)) {
            search_one_support(support, i, is_cs_and_coprime);
        }
    }
}


void fracessa::search_one_support(const bitset64& support, size_t support_size, bool is_cs_and_coprime)
{
    // `false` from the unsafe filter discards this support. `true` does not
    // accept anything; it only says that the exact test must make the decision.
    if (!conf_exact_) 
        if (!find_candidate_dbl(support, support_size))
            return;

    // This exact solve proves the Nash inequalities and fills candidate_.
    if (!find_candidate_frc(support, support_size))
            return;

    candidate_.support_size = support_size;
    candidate_.support = support;
    // IDs follow successful exact candidates in search order. Rejected support
    // masks never consume an ID.
    candidate_.candidate_id++;
    
    if (conf_with_log_)
        logger_->info("Found candidate! Check stability:");

    // The equilibrium condition is now proved; classify the stronger stability condition.
    check_stability();

    if (candidate_.is_ess)
        ess_count_++;

    // Rotated copies refer back to this solved representative.
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
    // Bomze's support criterion rules out an ESS on any strict superset of this
    // exact equilibrium support. Stability of the present candidate is irrelevant
    // to that pruning result.
    supports_.remove_supersets(candidate_.support, support_size);

    if (is_cs_and_coprime) {
        // Circular symmetry gives the remaining n-1 candidates without another
        // solve. Rotate the strategy and both support sets in the same direction.
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
