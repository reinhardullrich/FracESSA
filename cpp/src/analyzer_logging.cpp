#include <fracessa/fracessa.hpp>
#include <fracessa/fraction.hpp>

#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <sstream>
#include <string>
#include <string_view>

namespace fracessa {

namespace {

constexpr const char* search_method_name(search_method method) noexcept
{
    switch (method) {
        case search_method::fast: return "fast";
        case search_method::safe: return "safe";
    }
    return "unknown";
}

std::string rational_matrix_to_pretty_string(const numeric::matrix_int& matrix, numeric::integer::const_reference denominator)
{
    std::ostringstream stream;
    numeric::fraction value;
    for (size_t row = 0; row < matrix.rows(); ++row) {
        stream << "  " << row << ": [";
        for (size_t column = 0; column < matrix.cols(); ++column) {
            if (column != 0) stream << ", ";
            value.set_ratio(matrix(row, column), denominator);
            stream << value;
        }
        stream << ']';
        if (row + 1 != matrix.rows()) stream << '\n';
    }
    return stream.str();
}

} // namespace

void analyzer::start_log(search_method requested_method, bool is_cs, std::int64_t matrix_id)
{
    auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>("log/fracessa.log", 20 * 1024 * 1024, 5);
    logger_ = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
    logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
    logger_->set_level(spdlog::level::info);

    const std::string matrix_label = matrix_id >= 0 ? std::to_string(matrix_id) : "none";
    const std::string_view fallback = search::safe_fallback_name(safe_fallback_);
    constexpr std::string_view separator = "==============================================================================";
    logger_->info("");
    logger_->info("{}", separator);
    logger_->info("run started: matrix_id={} requested={} effective={} dimension={} circular={} fallback={}",
                  matrix_label, search_method_name(requested_method), search_method_name(method_), dimension_, is_cs,
                  fallback.empty() ? std::string_view{"none"} : fallback);
    logger_->info("{}", separator);
    logger_->info("game matrix ({} x {}):\n{}", game_.matrix.rows(), game_.matrix.cols(),
                  rational_matrix_to_pretty_string(game_.matrix, game_.denominator));
}

void analyzer::log_support_size(size_t support_size)
{
    if (!logger_ || support_size == last_logged_support_size_) return;
    logger_->info("searching support size: {}", support_size);
    last_logged_support_size_ = support_size;
}

void analyzer::log_candidate()
{
    if (!logger_) return;

    support_context_.copy(outside_best_replies_, candidate_.extended_support);
    support_context_.subtract(outside_best_replies_, candidate_.support);
    const size_t reference = support_context_.first(candidate_.support);

    logger_->info("solved candidate representative\n"
                  "  support:              {}\n"
                  "  extended:             {}\n"
                  "  reference:            {}\n"
                  "  outside best replies: {}",
                  support_context_.to_index_set(candidate_.support), support_context_.to_index_set(candidate_.extended_support), reference,
                  support_context_.to_index_set(outside_best_replies_));
}

void analyzer::log_reduced_b(const support::Support& outside_best_replies)
{
    if (!logger_) return;
    logger_->info("scaled reduced B for outside best replies {} ({} x {}):\n{}", support_context_.to_index_set(outside_best_replies),
                  scaled_reduced_b_.rows(), scaled_reduced_b_.cols(), scaled_reduced_b_.to_pretty_string());
}

void analyzer::finish_log()
{
    if (!logger_) return;
    logger_->info("run finished: weighted candidates={} weighted ESS={}", candidate_count_, ess_count_);
    logger_->flush();
}

} // namespace fracessa
