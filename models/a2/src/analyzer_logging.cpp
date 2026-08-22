#include <fracessa/fracessa.hpp>
#include <fracessa/bitset.hpp>
#include <fracessa/bitset_multiword.hpp>
#include <fracessa/fraction.hpp>

#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>

namespace fracessa {

namespace {

constexpr const char* search_method_name(search_method method) noexcept
{
    switch (method) {
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

template<class SupportMask>
void basic_analyzer<SupportMask>::start_log(search_method requested_method, bool is_cs, std::int64_t matrix_id)
{
    auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>("log/fracessa.log", 20 * 1024 * 1024, 5);
    logger_ = std::make_shared<spdlog::logger>("fracessa", rotating_sink);
    logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
    logger_->set_level(spdlog::level::info);

    const std::string matrix_label = matrix_id >= 0 ? std::to_string(matrix_id) : "none";
    constexpr std::string_view separator = "==============================================================================";
    logger_->info("");
    logger_->info("{}", separator);
    logger_->info("run started: matrix_id={} method={} dimension={} circular={}",
                  matrix_label, search_method_name(requested_method), dimension_, is_cs);
    logger_->info("{}", separator);
    logger_->info("game matrix ({} x {}):\n{}", game_.matrix.rows(), game_.matrix.cols(),
                  rational_matrix_to_pretty_string(game_.matrix, game_.denominator));
}

template<class SupportMask>
void basic_analyzer<SupportMask>::log_support_size(size_t support_size)
{
    if (!logger_ || support_size == last_logged_support_size_) return;
    logger_->info("searching support size: {}", support_size);
    last_logged_support_size_ = support_size;
}

template<class SupportMask>
void basic_analyzer<SupportMask>::log_candidate()
{
    if (!logger_) return;

    SupportMask outside_best_replies = candidate_.extended_support;
    size_t reference;
    if constexpr (std::is_same_v<SupportMask, support::bitset>) {
        outside_best_replies = support::subtract(outside_best_replies, candidate_.support);
        reference = support::find_pos_first_set_bit(candidate_.support);
    } else {
        outside_best_replies.subtract(candidate_.support);
        reference = candidate_.support.find_pos_first_set_bit();
    }

    logger_->info("solved candidate representative\n"
                  "  support:              {}\n"
                  "  extended:             {}\n"
                  "  reference:            {}\n"
                  "  outside best replies: {}",
                  support::to_index_set(candidate_.support), support::to_index_set(candidate_.extended_support), reference,
                  support::to_index_set(outside_best_replies));
}

template<class SupportMask>
void basic_analyzer<SupportMask>::log_reduced_b(const SupportMask& outside_best_replies)
{
    if (!logger_) return;
    logger_->info("scaled reduced B for outside best replies {} ({} x {}):\n{}", support::to_index_set(outside_best_replies),
                  scaled_reduced_b_.rows(), scaled_reduced_b_.cols(), scaled_reduced_b_.to_pretty_string());
}

template<class SupportMask>
void basic_analyzer<SupportMask>::finish_log()
{
    if (!logger_) return;
    logger_->info("run finished: weighted candidates={} weighted ESS={}", candidate_count_, ess_count_);
    logger_->flush();
}

template void basic_analyzer<support::bitset>::start_log(search_method, bool, std::int64_t);
template void basic_analyzer<support::bitset_multiword>::start_log(search_method, bool, std::int64_t);
template void basic_analyzer<support::bitset>::log_support_size(size_t);
template void basic_analyzer<support::bitset_multiword>::log_support_size(size_t);
template void basic_analyzer<support::bitset>::log_candidate();
template void basic_analyzer<support::bitset_multiword>::log_candidate();
template void basic_analyzer<support::bitset>::log_reduced_b(const support::bitset&);
template void basic_analyzer<support::bitset_multiword>::log_reduced_b(const support::bitset_multiword&);
template void basic_analyzer<support::bitset>::finish_log();
template void basic_analyzer<support::bitset_multiword>::finish_log();

} // namespace fracessa
