#include <fracessa/fracessa.hpp>
#include <fracessa/bitset64.hpp>
#include <linalg/matrix_fraction.hpp>

#include <utility>

namespace {

void rotate_candidate(candidate& value, size_t dimension) {
    // Apply the same one-position rotation to both support masks and x.
    value.support = bs64::rot_one_right(value.support, dimension);
    value.extended_support =
        bs64::rot_one_right(value.extended_support, dimension);

    fraction first = value.vector(0, 0);
    for (size_t i = 0; i + 1 < dimension; ++i)
        value.vector(i, 0) = value.vector(i + 1, 0);
    value.vector(dimension - 1, 0) = first;
}

void reflect_candidate(candidate& value, size_t dimension) {
    // Mirror strategy i to n-1-i in the masks and in x.
    value.support = reflect_support(value.support, dimension);
    value.extended_support = reflect_support(value.extended_support, dimension);
    for (size_t i = 0; i < dimension / 2; ++i)
        std::swap(value.vector(i, 0), value.vector(dimension - 1 - i, 0));
}

} // namespace

fracessa::fracessa(
    const linalg::matrix_frc& matrix, bool is_cs, bool with_candidates,
    bool exact, bool full_support, bool with_log, int matrix_id)
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
            "log/fracessa.log", 20 * 1024 * 1024, 5);
        logger_ = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
        logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
        logger_->set_level(spdlog::level::info);

        logger_->info("");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");
        logger_->info("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#");
        logger_->info("*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*");
        if (matrix_id >= 0)
            logger_->info("matrix_id={}", matrix_id);

        logger_->info("n={}", dimension_);
        logger_->info(
            "game matrix:\n{}",
            matrix_server_.get_game_matrix_frc().to_log_string());
    }

    if (!conf_exact_ && !matrix_server_.initialize_unsafe_filter())
        conf_exact_ = true;

    if (conf_full_support_) {
        const bitset64 full_support_mask = bs64::set_all_n_bits(dimension_);
        search_one_support(full_support_mask, dimension_);
        if (ess_count_ > 0)
            return;
    }

    supports_.initialize();

    const size_t max_support_size =
        conf_full_support_ ? dimension_ - 1 : dimension_;
    for (size_t support_size = 1;
         support_size <= max_support_size; ++support_size) {
        if (conf_with_log_)
            logger_->info("Searching support size {}:", support_size);

        for (const bitset64 support : supports_.get_supports(support_size))
            search_one_support(support, support_size);
    }
}

void fracessa::search_one_support(
    const bitset64& support, size_t support_size, bool)
{
    if (!conf_exact_ && !find_candidate_dbl(support, support_size))
        return;

    if (!find_candidate_frc(support, support_size))
        return;

    candidate_.support_size = support_size;
    candidate_.support = support;

    if (conf_with_log_)
        logger_->info("Found candidate! Check stability:");

    check_stability();

    /*
     * This exact equilibrium support excludes all its strict supersets. In a
     * circular game the same is true after every rotation and reflection, so
     * remove_supersets() records the complete family as future pruning rules.
     * It also returns the family's shape, which lets us count and output all
     * equivalent candidates without solving their linear systems separately.
     */
    const SupportOrbit orbit =
        supports_.remove_supersets(candidate_.support, support_size);
    const size_t orbit_size = orbit.member_count();

    if (candidate_.is_ess)
        ess_count_ += orbit_size;

    const size_t first_id = candidate_.candidate_id + 1;
    candidate_.shift_reference = orbit_size > 1 ? first_id : 0;

    if (!conf_with_candidates_ && !conf_with_log_) {
        // No rows are needed, but IDs still count every equivalent candidate.
        candidate_.candidate_id += orbit_size;
        return;
    }

    // Save the solved orientation before the first rotation; the reflected
    // family, when distinct, must start from this original candidate.
    candidate reflected_base;
    if (orbit.has_reflected_class)
        reflected_base = candidate_;

    size_t next_id = first_id;
    const auto emit = [&]() {
        candidate_.candidate_id = next_id++;
        if (conf_with_candidates_)
            candidates_.push_back(candidate_);
        if (conf_with_log_) {
            logger_->info("{}", candidate::header());
            logger_->info("{}", candidate_.to_string());
        }
    };

    for (size_t rotation = 0;
         rotation < orbit.rotation_count; ++rotation) {
        if (rotation != 0)
            rotate_candidate(candidate_, dimension_);
        emit();
    }

    if (orbit.has_reflected_class) {
        candidate_ = std::move(reflected_base);
        reflect_candidate(candidate_, dimension_);
        for (size_t rotation = 0;
             rotation < orbit.rotation_count; ++rotation) {
            if (rotation != 0)
                rotate_candidate(candidate_, dimension_);
            emit();
        }
    }

    candidate_.candidate_id = next_id - 1;
}
