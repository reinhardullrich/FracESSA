#pragma once

#include <cstddef>
#include <cstdint>

#include <fracessa/bitset64.hpp>
#include <linalg/fraction.hpp>
#include <linalg/matrix_double.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

// The arithmetic helpers are declared here so the focused tests exercise the
// same strictly compiled implementation as production.

struct rational_bounds {
    double midpoint;
    double radius;
    double magnitude;
};

struct solution_error_bound {
    double q;
    double error;
};

class find_candidate_verified {
public:
    explicit find_candidate_verified(const linalg::matrix_frc& game_matrix);

    // False means the support is proven not to be a candidate. True means the
    // result is inconclusive and exact arithmetic must decide.
    bool find(const bitset64& support, size_t support_size);

private:
    void initialize_bounds();

    const linalg::matrix_frc& game_frc_;
    size_t dimension_;
    bool bounds_initialized_;
    bool bounds_available_;
    linalg::matrix_dbl game_midpoint_;
    linalg::matrix_dbl game_radius_;
    linalg::matrix_dbl game_magnitude_;
    linalg::matrix_dbl compact_lu_;
};

const char* unavailable_reason() noexcept;
bool round_down(double value, double& result) noexcept;
bool round_up(double value, double& result) noexcept;
bool get_rational_bounds(const fraction& value, rational_bounds& result) noexcept;

bool factor_lu(linalg::matrix_dbl& matrix, size_t dimension, uint8_t* permutation) noexcept;
bool triangular_solve_in_place(const linalg::matrix_dbl& matrix, size_t dimension, double* solution) noexcept;
bool absolute_lu_apply(
    const linalg::matrix_dbl& matrix,
    size_t dimension,
    const double* input,
    double* output
) noexcept;

bool prove_solution_error(const linalg::matrix_dbl& matrix, size_t dimension, const uint8_t* permutation,
                          const double* input_row_bound, const double* residual_lower, const double* residual_upper,
                          solution_error_bound& result) noexcept;

bool residual_enclosure(const linalg::matrix_dbl& game_midpoint, const linalg::matrix_dbl& game_radius,
                        const uint8_t* support_indices, size_t support_size, const uint8_t* permutation,
                        const double* solution, double* input_row_bound, double* residual_lower,
                        double* residual_upper) noexcept;

bool outside_gain_lower_bound(const linalg::matrix_dbl& game_midpoint, const linalg::matrix_dbl& game_radius,
                              const linalg::matrix_dbl& game_magnitude, size_t outside_index,
                              const uint8_t* support_indices, size_t support_size, const double* solution,
                              double error, double& gain_lower) noexcept;

} // namespace candidate_search
