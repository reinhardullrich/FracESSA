#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <rational_linalg/matrix.hpp>
#include <fracessa/bitset64.hpp>
#include <string>
#include <sstream>
#include <iomanip>

class candidate
{
    public:
        size_t candidate_id = 0;
        rational_linalg::Matrix<fraction> vector;  // Column vector: Matrix<fraction>(n, 1) - always use arbitrary precision
        bitset64 support;
        size_t support_size = 0;
        bitset64 extended_support;
        size_t extended_support_size;
        size_t shift_reference;
        bool is_ess;
        std::string stability;
        fraction payoff;  // Always use arbitrary precision
        double payoff_double;

        std::string to_string() const
        {
            std::ostringstream oss;
            oss << candidate_id << ";";
            for (size_t i = 0; i < vector.rows(); i++) {
                oss << vector(i, 0);  // Uses optimized operator<<, no string allocation
                if (i < vector.rows() - 1) {
                    oss << ",";
                }
            }
            oss << ";" << bs64::to_string(support) << ";"
                << support_size << ";"
                << bs64::to_string(extended_support) << ";"
                << extended_support_size << ";"
                << shift_reference << ";"
                << is_ess << ";"
                << stability << ";"
                << payoff << ";"  // Uses optimized operator<<, no string allocation
                << std::fixed << std::setprecision(6) << payoff_double;
            return oss.str();
        }

        static std::string header()
        {
            return "candidate_id;vector;support;support_size;extended_support;extended_support_size;shift_reference;is_ess;stability;payoff;payoff_double;";
        }
};

#endif // CANDIDATE_H
