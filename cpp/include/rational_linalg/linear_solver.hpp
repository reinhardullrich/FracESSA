#ifndef RATIONAL_LINALG_LINEAR_SOLVER_HPP
#define RATIONAL_LINALG_LINEAR_SOLVER_HPP

#include <rational_linalg/matrix.hpp>
#include <rational_linalg/constants.hpp>
#include <cmath>
#include <limits>

namespace rational_linalg {

/*
 * BareissGauss (fraction-free Gaussian elimination) solver for Matrix<T>
 *
 * Speed:
 *   - Complexity: O(n^3), same as standard Gaussian elimination.
 *   - Much slower in practice with fractional or big integer types
 *     because intermediate numbers grow rapidly.
 *   - For floating-point, Bareiss is slower than standard methods
 *     because it avoids floating-point division tricks and uses exact arithmetic.
 *
 * Stability / Accuracy:
 *   - Perfect for fractional or symbolic types; no rounding errors.
 *   - Detects singular matrices exactly (pivot = 0).
 *   - Great for ill-conditioned matrices in exact arithmetic because
 *     no division by small numbers until necessary.
 *   - Cannot suffer from numerical instability in exact arithmetic.
 *
 * Use case:
 *   - Exact solutions, guaranteed singularity detection,
 *     symbolic or fractional arithmetic.
 */

template<typename T>
class BareissGauss {
public:
    using Scalar = T;
    using VectorType = Matrix<T>;

    explicit BareissGauss(const Matrix<T>& Ab) : M(Ab) {}

    bool solve(Matrix<T>& x) {
        const size_t n = M.rows();

        T divPrev = rational_linalg::one<T>();

        // Bareiss elimination with partial pivoting
        for (size_t k = 0; k < n - 1; ++k) {
            // Partial pivoting: find row with maximum absolute value in column k
            size_t max_row = k;
            T max_val = M(k, k);
            if (max_val < rational_linalg::zero<T>()) max_val = -max_val;
            
            for (size_t i = k + 1; i < n; ++i) {
                T val = M(i, k);
                if (val < rational_linalg::zero<T>()) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            // Swap rows if necessary
            if (max_row != k) {
                M.swap_rows(k, max_row);
            }
            
            const T pivot = M(k, k);
            if (pivot == rational_linalg::zero<T>()) {
                return false;
            }

            // Bareiss elimination step
            for (size_t i = k + 1; i < n; ++i) {
                for (size_t j = k + 1; j <= n; ++j) {
                    M(i, j) = (M(i, j) * pivot - M(i, k) * M(k, j)) / divPrev;
                }
                M(i, k) = rational_linalg::zero<T>(); // Zero out eliminated element (not updated in inner loop)
            }

            divPrev = pivot;
        }

        if (M(n - 1, n - 1) == rational_linalg::zero<T>()) return false;

        // Back substitution
        x = Matrix<T>(n, 1);
        for (size_t i = n; i-- > 0; ) {
            T sum = M(i, n);
            for (size_t j = i + 1; j < n; ++j) {
                sum -= M(i, j) * x(j, 0);
            }

            const T pivot = M(i, i);

            if (pivot == rational_linalg::zero<T>()) 
                return false;            

            T temp_x = sum / pivot;
            if (temp_x <= rational_linalg::zero<T>())
                return false; //early exit if sum is not >zero! then it cannot be a solution (since x on the simplex is always positive!)
            else
                x(i, 0) = temp_x;
        }

        return true;
    }

private:
    Matrix<T> M;
};

/*
 * GaussDouble - Optimized standard Gaussian elimination for Matrix<double>
 *
 * Speed:
 *   - Complexity: O(n^3), same as BareissGauss
 *   - Much faster than BareissGauss for double because:
 *     - Uses direct division (no fraction-free arithmetic overhead)
 *     - Optimized for floating-point operations
 *     - Better cache locality
 *
 * Stability / Accuracy:
 *   - Uses partial pivoting for numerical stability
 *   - Tolerance-based singularity detection
 *   - Standard floating-point precision limitations apply
 *
 * Use case:
 *   - Fast solving of double-precision linear systems
 *   - When exact arithmetic is not required
 */


// class GaussDouble {
// public:
//     using Scalar = double;
//     using VectorType = Matrix<double>;

//     explicit GaussDouble(const Matrix<double>& Ab) : M(Ab) {}

//     bool solve(Matrix<double>& x) {
//         const size_t n = M.rows();

//         const double errorbound = 1e-5 * n; // huge margin, false positives eliminated by fractional check

//         const double epsilon = 1e-15; //std::numeric_limits<double>::epsilon(); // * n; keep epsilon super small to avoid false negative. false positive will be detected by fractional afterwards!!!

//         // Standard Gaussian elimination with partial pivoting
//         for (size_t k = 0; k < n - 1; ++k) {
//             // Partial pivoting: find row with maximum absolute value in column k
//             size_t max_row = k;
//             double max_val = std::abs(M(k, k));
            
//             for (size_t i = k + 1; i < n; ++i) {
//                 double val = std::abs(M(i, k));
//                 if (val > max_val) {
//                     max_val = val;
//                     max_row = i;
//                 }
//             }
            
//             // Check for singularity
//             if (max_val < epsilon) {
//                 return false;
//             }
            
//             // Swap rows if necessary
//             if (max_row != k) {
//                 M.swap_rows(k, max_row);
//             }
            
//             const double pivot = M(k, k);

//             // Standard elimination step (direct division, no fraction-free)
//             for (size_t i = k + 1; i < n; ++i) {
//                 const double factor = M(i, k) / pivot;
//                 for (size_t j = k + 1; j <= n; ++j) {
//                     M(i, j) -= factor * M(k, j);
//                 }
//                 M(i, k) = 0.0; // Zero out eliminated element (not updated in inner loop)
//             }
//         }

//         // Check last pivot
//         if (std::abs(M(n - 1, n - 1)) < epsilon) {
//             return false;
//         }

//         // Back substitution
//         x = Matrix<double>(n, 1);
//         for (size_t i = n; i-- > 0; ) {
//             double sum = M(i, n);
//             for (size_t j = i + 1; j < n; ++j) {
//                 sum -= M(i, j) * x(j, 0);
//             }

//             const double pivot = M(i, i);
//             if (std::abs(pivot) < epsilon) {
//                 return false;
//             }

//             double temp_x = sum / pivot;
//             // Check for not >zero in solution component with huge margin
//             // False positives will be eliminated by fractional checks later
//             if (temp_x < -errorbound) 
//                 return false;
//             else
//                 x(i, 0) = temp_x;
//         }

//         return true;
//     }

// private:
//     Matrix<double> M;
// };

/*
 * GaussDouble - STRICTLY OPTIMIZED
 * 
 * 1. NO HEAP ALLOCATIONS (Fixes the slowdown).
 * 2. Pointer Hoisting: Calculates row addresses once per row.
 * 3. Physical Swapping: Faster than virtual pivoting for small N (fits in L1 cache).
 */
 class GaussDouble {
    public:
        using Scalar = double;
        using VectorType = Matrix<double>;
    
        explicit GaussDouble(Matrix<double>& Ab) : M(Ab) {}
    
        bool solve(Matrix<double>& x) {
            const size_t n = M.rows();
            const size_t cols = M.cols(); // n+1 usually
            double* data = M.data();
    
            const double errorbound = 1e-5 * n; 
            const double epsilon = 1e-15; 
    
            // Forward elimination
            for (size_t k = 0; k < n - 1; ++k) {
                
                // 1. Partial Pivoting
                size_t max_row = k;
                double max_val = std::abs(data[k * cols + k]);
                
                for (size_t i = k + 1; i < n; ++i) {
                    double val = std::abs(data[i * cols + k]);
                    if (val > max_val) {
                        max_val = val;
                        max_row = i;
                    }
                }
                
                if (max_val < epsilon) return false;
                
                // Physical Swap (Fast for small N)
                if (max_row != k) {
                    // Swap rows manually using pointers to avoid index math overhead
                    double* row_k = data + k * cols;
                    double* row_max = data + max_row * cols;
                    // Swap the relevant part of the row (k to end is enough, but cols is safe)
                    // We swap the whole row to be safe with the layout
                    for(size_t j=0; j < cols; ++j) {
                         std::swap(row_k[j], row_max[j]);
                    }
                }
                
                // Hoist pivot pointer
                double* pivot_row = data + k * cols;
                const double pivot = pivot_row[k];
    
                // 2. Elimination
                for (size_t i = k + 1; i < n; ++i) {
                    double* curr_row = data + i * cols;
                    
                    const double factor = curr_row[k] / pivot;
                    curr_row[k] = 0.0; 
                    
                    // Vectorizable inner loop
                    for (size_t j = k + 1; j <= n; ++j) { // <= n includes the b vector
                        curr_row[j] -= factor * pivot_row[j];
                    }
                }
            }
    
            // Check last pivot
            if (std::abs(data[(n - 1) * cols + (n - 1)]) < epsilon) {
                return false;
            }
    
            // 3. Back substitution
            x = Matrix<double>(n, 1);
            for (size_t i = n; i-- > 0; ) {
                const double* row_ptr = data + i * cols;
                double sum = row_ptr[n]; // b vector
                
                for (size_t j = i + 1; j < n; ++j) {
                    sum -= row_ptr[j] * x(j, 0);
                }
    
                const double pivot = row_ptr[i];
                if (std::abs(pivot) < epsilon) return false;
    
                double temp_x = sum / pivot;
                if (temp_x < -errorbound) return false;
                
                x(i, 0) = temp_x;
            }
    
            return true;
        }
    
    private:
        Matrix<double>& M;
    };


    
/*
 * GaussRational<T> - Generic standard Gaussian elimination for Matrix<T>
 *
 * Speed:
 *   - Complexity: O(n^3), same as BareissGauss
 *   - Faster than BareissGauss for types that support efficient division:
 *     - Uses direct division (no fraction-free arithmetic overhead)
 *     - Better for fractional types when intermediate values don't grow too large
 *     - More efficient for floating-point types
 *
 * Stability / Accuracy:
 *   - Uses partial pivoting for numerical stability
 *   - Exact zero detection for singularity (no tolerance needed for exact types)
 *   - For floating-point types, may suffer from numerical instability
 *     compared to BareissGauss for ill-conditioned matrices
 *
 * Use case:
 *   - Fast solving of linear systems when exact arithmetic is available
 *   - Alternative to BareissGauss when division is efficient
 *   - For fractional types when matrix is well-conditioned
 */
// template<typename T>
// class GaussRational {
// public:
//     using Scalar = T;
//     using VectorType = Matrix<T>;

//     explicit GaussRational(const Matrix<T>& Ab) : M(Ab) {}

//     bool solve(Matrix<T>& x) {
//         const size_t n = M.rows();

//         // Standard Gaussian elimination with partial pivoting
//         for (size_t k = 0; k < n - 1; ++k) {
//             // Partial pivoting: find row with maximum absolute value in column k
//             size_t max_row = k;
//             T max_val = M(k, k);
//             if (max_val < T(0)) max_val = -max_val;
            
//             for (size_t i = k + 1; i < n; ++i) {
//                 T val = M(i, k);
//                 if (val < T(0)) val = -val;
//                 if (val > max_val) {
//                     max_val = val;
//                     max_row = i;
//                 }
//             }
            
//             // Check for singularity (exact zero check for generic types)
//             if (max_val == T(0)) {
//                 return false;
//             }
            
//             // Swap rows if necessary
//             if (max_row != k) {
//                 M.swap_rows(k, max_row);
//             }
            
//             const T pivot = M(k, k);

//             // Standard elimination step (direct division, no fraction-free)
//             for (size_t i = k + 1; i < n; ++i) {
//                 const T factor = M(i, k) / pivot;
//                 for (size_t j = k + 1; j <= n; ++j) {
//                     M(i, j) -= factor * M(k, j);
//                 }
//                 M(i, k) = T(0); // Zero out eliminated element (not updated in inner loop)
//             }
//         }

//         // Check last pivot
//         if (M(n - 1, n - 1) == T(0)) {
//             return false;
//         }

//         // Back substitution
//         x = Matrix<T>(n, 1);
//         for (size_t i = n; i-- > 0; ) {
//             T sum = M(i, n);
//             for (size_t j = i + 1; j < n; ++j) {
//                 sum -= M(i, j) * x(j, 0);
//             }

//             const T pivot = M(i, i);

//             if (pivot == T(0)) 
//                 return false;

//             T temp_x = sum / pivot;
//             if (temp_x <= T(0))
//                 return false; //early exit if sum is not >zero! 
//                 //then it cannot be a solution: full support for this subgame <=> x is in the interior of the (sub)simplex <=> all x_i are positive!
//             else
//                 x(i, 0) = temp_x;               
//         }

//         return true;
//     }

// private:
//     Matrix<T> M;
// };

/*
 * GaussRational<T> - Optimized Standard Gaussian elimination.
 * 
 * OPTIMIZATIONS:
 * 1. Uses fmpq_submul (FLINT) for direct subtraction-multiplication without temporaries.
 * 2. Uses Virtual Pivoting to avoid expensive row swaps in memory.
 */
 template<typename T>
 class GaussRational {
 public:
     using Scalar = T;
     using VectorType = Matrix<T>;
 
     // Constructor accepts non-const ref to modify directly (copy is made by caller/MatrixServer usually, or we assume local scratchpad)
     explicit GaussRational(Matrix<T>& Ab) : M(Ab) {}
 
     bool solve(Matrix<T>& x) {
         const size_t n = M.rows();
 
         // Virtual pivoting map: p[i] stores the physical row index for logical row i
         std::vector<size_t> p(n);
         for (size_t i = 0; i < n; ++i) p[i] = i;
 
         // Forward elimination
         for (size_t k = 0; k < n - 1; ++k) {
             // Partial pivoting
             size_t max_row_logical = k;
             
             // For exact types, finding ANY non-zero is sufficient for solvability, 
             // though max abs value keeps coefficient growth slightly cleaner.
             // Using a simple non-zero check is fastest if we trust FLINT's GCD.
             // Let's stick to simple "first non-zero" or "current if non-zero" to avoid abs() overhead, 
             // or implement full abs comparison if needed. 
             // Here, we check if current pivot is zero. If so, search for a non-zero.
             
             fmpq* max_val_ptr = M(p[k], k).ptr();
             
             if (fmpq_is_zero(max_val_ptr)) {
                 bool found = false;
                 for (size_t i = k + 1; i < n; ++i) {
                     fmpq* val_ptr = M(p[i], k).ptr();
                     if (!fmpq_is_zero(val_ptr)) {
                         max_row_logical = i;
                         max_val_ptr = val_ptr;
                         found = true;
                         break;
                     }
                 }
                 if (!found) return false; // Singular column
             }
             
             // Swap virtual rows
             if (max_row_logical != k) {
                 std::swap(p[k], p[max_row_logical]);
             }
             
             const size_t pivot_row_idx = p[k];
             fmpq* pivot_ptr = M(pivot_row_idx, k).ptr(); // Corrected type: fmpq*
 
             // Elimination
             for (size_t i = k + 1; i < n; ++i) {
                 size_t curr_row_idx = p[i];
                 fmpq* target_ptr = M(curr_row_idx, k).ptr(); // Corrected type: fmpq*
 
                 if (fmpq_is_zero(target_ptr)) continue;
 
                 // factor = M(i, k) / pivot
                 // We use a temporary fmpq_t. It acts as fmpq* when passed to functions.
                 fmpq_t factor;
                 fmpq_init(factor);
                 fmpq_div(factor, target_ptr, pivot_ptr);
 
                 // Row_i -= factor * Row_k
                 // Optimize: Start from k+1 (cols 0..k become 0)
                 for (size_t j = k + 1; j <= n; ++j) {
                     // M(i, j) = M(i, j) - factor * M(k, j)
                     // fmpq_submul(rop, op1, op2) -> rop = rop - op1 * op2
                     fmpq_submul(M(curr_row_idx, j).ptr(), factor, M(pivot_row_idx, j).ptr());
                 }
                 
                 // Explicitly zero out the eliminated element
                 fmpq_zero(target_ptr);
                 fmpq_clear(factor);
             }
         }
 
         // Check last pivot
         if (M(p[n - 1], n - 1).is_zero()) {
             return false;
         }
 
         // Back substitution
         x = Matrix<T>(n, 1);
         for (size_t i = n; i-- > 0; ) {
             size_t row_idx = p[i];
             
             // sum = b[i] (which is in column n)
             // Use temp to accumulate to safely use fmpq_submul
             fmpq_t sum;
             fmpq_init(sum);
             fmpq_set(sum, M(row_idx, n).ptr());
 
             for (size_t j = i + 1; j < n; ++j) {
                 fmpq_submul(sum, M(row_idx, j).ptr(), x(j, 0).ptr());
             }
 
             // x[i] = sum / M[i, i]
             fmpq* pivot_ptr = M(row_idx, i).ptr();
             
             // Division
             fmpq_div(x(i, 0).ptr(), sum, pivot_ptr);
             
             // Early exit: x must be > 0
             // fmpq_sgn returns positive (>0), zero (0), negative (<0)
             if (fmpq_sgn(x(i, 0).ptr()) <= 0) {
                 fmpq_clear(sum);
                 return false; 
             }
 
             fmpq_clear(sum);
         }
 
         return true;
     }
 
 private:
     Matrix<T>& M;
 };


// Helper struct for automatic solver selection
// Automatically selects GaussDouble for double, GaussRational<T> for other types
template<typename T>
struct SolverSelector {
    using type = GaussRational<T>;
    //using type = BareissGauss<T>; //leave this here!!! for comparison!!!
};

template<>
struct SolverSelector<double> {
    using type = GaussDouble;
};

// Type alias for automatic solver selection
template<typename T>
using LinearSolver = typename SolverSelector<T>::type;

} // namespace rational_linalg

#endif // RATIONAL_LINALG_LINEAR_SOLVER_HPP

