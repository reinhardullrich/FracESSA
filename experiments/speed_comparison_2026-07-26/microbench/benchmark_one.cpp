#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fracessa/fracessa.hpp>
#include <fracessa/matrix_parser.hpp>

namespace {

using clock_type = std::chrono::steady_clock;

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    if ((values.size() & 1U) != 0U) {
        return values[middle];
    }
    return (values[middle - 1] + values[middle]) / 2.0;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " MATRIX [TARGET_SECONDS]\n";
        return 2;
    }

    const double target_seconds = argc == 3 ? std::stod(argv[2]) : 3.0;
    if (!(target_seconds > 0.0) || !std::isfinite(target_seconds)) {
        throw std::invalid_argument("TARGET_SECONDS must be finite and positive");
    }

    linalg::matrix_frc matrix;
    bool is_cs = false;
    if (!matrix_parser::parse_matrix_string(argv[1], matrix, is_cs)) {
        std::cerr << "Invalid matrix\n";
        return 2;
    }

    std::uint64_t checksum = 0;
    size_t last_ess_count = 0;
    size_t last_candidate_count = 0;

    auto execute_once = [&]() {
        ::fracessa analyzer(matrix, is_cs, true, false, false, false, -1);
        last_ess_count = analyzer.ess_count_;
        last_candidate_count = analyzer.candidates_.size();
        checksum = checksum * 1315423911ULL + last_ess_count * 257ULL + last_candidate_count;
    };

    auto measure = [&](size_t repetitions) {
        const auto start = clock_type::now();
        for (size_t i = 0; i < repetitions; ++i) {
            execute_once();
        }
        const auto end = clock_type::now();
        return std::chrono::duration<double, std::nano>(end - start).count();
    };

    const double target_ns = target_seconds * 1.0e9;
    const double pilot_ns = std::max(1.0, measure(1));
    size_t warmup_repetitions = 0;
    size_t calibration_repetitions = 0;
    double post_warm_pilot_ns = pilot_ns;
    double calibrated_ns = pilot_ns;
    size_t measured_repetitions = 1;
    std::vector<double> per_run_samples;
    double measured_total_ns = pilot_ns;

    if (pilot_ns < target_ns) {
        constexpr double warmup_target_ns = 50.0e6;
        if (pilot_ns < warmup_target_ns) {
            warmup_repetitions = static_cast<size_t>(std::ceil(warmup_target_ns / pilot_ns));
            for (size_t i = 0; i < warmup_repetitions; ++i) {
                execute_once();
            }
        }

        const double calibration_target_ns = std::min(20.0e6, target_ns * 0.05);
        post_warm_pilot_ns = std::max(1.0, measure(1));
        calibration_repetitions = static_cast<size_t>(std::ceil(calibration_target_ns / post_warm_pilot_ns));
        calibrated_ns = measure(calibration_repetitions) / static_cast<double>(calibration_repetitions);
        measured_repetitions = static_cast<size_t>(std::ceil(target_ns / calibrated_ns));
        const size_t batch_count = std::min<size_t>(7, measured_repetitions);
        const size_t base = measured_repetitions / batch_count;
        const size_t remainder = measured_repetitions % batch_count;

        per_run_samples.reserve(batch_count);
        measured_total_ns = 0.0;
        for (size_t batch = 0; batch < batch_count; ++batch) {
            const size_t batch_repetitions = base + (batch < remainder ? 1 : 0);
            const double batch_ns = measure(batch_repetitions);
            measured_total_ns += batch_ns;
            per_run_samples.push_back(batch_ns / static_cast<double>(batch_repetitions));
        }
    } else {
        per_run_samples.push_back(pilot_ns);
    }

    const double mean_ns = measured_total_ns / static_cast<double>(measured_repetitions);
    const auto [minimum, maximum] = std::minmax_element(per_run_samples.begin(), per_run_samples.end());

    std::cout << std::setprecision(17)
              << "repetitions=" << measured_repetitions << '\n'
              << "batches=" << per_run_samples.size() << '\n'
              << "warmup_repetitions=" << warmup_repetitions << '\n'
              << "calibration_repetitions=" << calibration_repetitions << '\n'
              << "pilot_ns=" << pilot_ns << '\n'
              << "post_warm_pilot_ns=" << post_warm_pilot_ns << '\n'
              << "calibrated_ns=" << calibrated_ns << '\n'
              << "median_ns=" << median(per_run_samples) << '\n'
              << "mean_ns=" << mean_ns << '\n'
              << "min_ns=" << *minimum << '\n'
              << "max_ns=" << *maximum << '\n'
              << "measured_seconds=" << measured_total_ns / 1.0e9 << '\n'
              << "ess_count=" << last_ess_count << '\n'
              << "candidate_count=" << last_candidate_count << '\n'
              << "checksum=" << checksum << '\n';
    return 0;
}
