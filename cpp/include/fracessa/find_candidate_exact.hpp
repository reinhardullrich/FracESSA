#pragma once

#include <cstddef>

#include <fracessa/bitset64.hpp>
#include <fracessa/candidate.hpp>
#include <linalg/matrix_fraction.hpp>

namespace candidate_search {

class find_candidate_exact {
public:
    explicit find_candidate_exact(const linalg::matrix_frc& game_matrix) noexcept;

    // True means an exact candidate was found and written to result. The dense vector is optional because stability does not use it.
    bool find(const bitset64& support, size_t support_size, candidate& result, bool materialize_vector);

    // Valid after find() returns true. The reduced Hessian is Z^T*A_S*Z, where the columns of Z span the zero-sum directions on the support.
    // Its negative definiteness is the support-only second-order ESS condition.
    bool reduced_hessian_is_negative_definite() const noexcept { return reduced_hessian_is_negative_definite_; }

private:
    const linalg::matrix_frc& game_frc_;
    size_t dimension_;
    // Reused exact workspaces for the reduced system and its reconstructed probability vector. Their dimensions change only with support size.
    linalg::matrix_frc reduced_system_;
    linalg::matrix_frc support_solution_;
    bool reduced_hessian_is_negative_definite_ = false;
};

} // namespace candidate_search
