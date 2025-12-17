#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <rational_linalg/matrix.hpp>
#include <fracessa/bitset64.hpp>
#include <string>

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

        std::string to_string()
        {
            std::string str;
            str.reserve(256);  // Pre-allocate reasonable size
            str += std::to_string(candidate_id) + ";";
            for (size_t i = 0; i < vector.rows(); i++) {
                // Convert fraction to string, normalize integers (X/1 -> X)
                str += vector(i, 0).to_string() + ",";
            }
            if (vector.rows() > 0)
                str.pop_back();
            str += ";" + bs64::to_string(support) + ";";
            str += std::to_string(support_size) + ";";
            str += bs64::to_string(extended_support) + ";";
            str += std::to_string(extended_support_size) + ";";
            str += std::to_string(shift_reference) + ";";
            str += std::to_string(is_ess) + ";";
            str += stability + ";";
            str += payoff.to_string() + ";";
            str += std::to_string(payoff_double);

            return str;
        }

        static std::string header()
        {
            return "candidate_id;vector;support;support_size;extended_support;extended_support_size;shift_reference;is_ess;stability;payoff;payoff_double;";
        }
};

#endif // CANDIDATE_H
