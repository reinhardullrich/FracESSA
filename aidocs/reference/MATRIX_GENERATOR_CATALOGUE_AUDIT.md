# Matrix Generator Catalogue Audit

## Purpose

This audit treats Anymatrix, TypedMatrices.jl, and Matrix Depot as one overlapping matrix-generator family. It imports a broad,
deterministic sample of matrices useful for FracESSA while enforcing the project contract:

1. The matrix is square and symmetric.
2. Every entry is an exact integer or rational number.
3. The dimension is between 1 and 63 inclusive.
4. An exact matrix already in the database is not inserted again.
5. A circulant matrix is stored as the strategically equivalent zero-diagonal game `A - dJ`; dimension one is the documented
   compact-format exception.
6. Catalog imports have null candidate and ESS baselines until they are analyzed explicitly.

## Audited Sources

| Catalogue | Revision | Repository | License in audited revision |
| --- | --- | --- | --- |
| Anymatrix | `fc8065712db8aaf01d6c09354b4824da7704747a` | https://github.com/north-numerical-computing/anymatrix | BSD-2-Clause |
| TypedMatrices.jl | `5a0b08f8e86fe218673fbc7c3c8636759b9e7b33` | https://github.com/TypedMatrices/TypedMatrices.jl | MIT |
| Matrix Depot | `7bda04519d0b0fc10249b31ce885231f14ffa7d1` | https://github.com/JuliaMatrices/MatrixDepot.jl | MIT |

The database stores a revision-pinned source URL and original generator identifier for each imported matrix. Overlapping formulas are
merged before sampling. In particular, Dingdong and RIS use the same exact formula and are represented as one family.

## Sampling Rule

The five selection bands are:

| Band | Persistent `size_class` |
| --- | --- |
| 1-8 | `small` |
| 9-16 | `medium` |
| 17-25 | `large` |
| 26-44 | `super_large` |
| 45-63 | `super_large` |

For each populated `property x band` stratum, candidates are ranked by SHA-256 using seed `20260802`, and at most three distinct exact
matrices are retained. A second deterministic pass adds a representative for any eligible family not reached by property sampling.
One matrix may satisfy several strata or families, so the selected total is much smaller than the sum of all stratum quotas.

The 34 classifications are:

`binary`, `circulant`, `diagonal`, `diagonally_dominant`, `eigenvalues`, `fixed_size`, `graph`, `hankel`, `ill_conditioned`,
`indefinite`, `infinitely_divisible`, `integer`, `inverse`, `involutory`, `m_matrix`, `nonnegative`, `normal`, `orthogonal`,
`parameter_dependent`, `positive`, `positive_definite`, `positive_semidefinite`, `random`, `rank_deficient`,
`regularisation_problem`, `sparse`, `stochastic`, `symmetric`, `toeplitz`, `totally_nonnegative`, `totally_positive`,
`transformation`, `tridiagonal`, and `unimodular`.

These tags combine source catalogue classifications with exact structural facts such as symmetry, integrality, sparsity, or circulant
form. They are selection metadata, not newly proved general theorems about each named generator family.

## Eligible Families

The audit found 51 eligible exact families:

`augment`, `beta`, `block_householder`, `cauchy`, `circulant`, `circulant_binomial`, `cross`, `deriv2`, `dingdong_ris`, `fiedler`,
`gcd`, `hankel`, `hilbert`, `inverse_hilbert`, `indefinite_power_sum`, `ipjfact_factorial`, `ipjfact_reciprocal`, `kms`, `lehmer`,
`minij`, `moler`, `pascal`, `pei`, `poisson`, `rosser`, `strakos`, `symmetric_stochastic`, `toeplitz`, `tridiagonal`, `wilkinson`,
`zielke_symmetric`, `dembo9`, `wilson`, `hexgrid`, `erdos_renyi`, `gilbert`, `geometric_random_graph`, `lock_and_key`,
`preferential_attachment`, `range_dependent_graph`, `small_world`, `sticky_graph`, `graph_laplacian`, `graph_path_length`,
`rewired_graph`, `shortcut_graph`, `bait_prey_sample`, `uniform_node_sample`, `kleinberg`, `benguela`, and `hadamard`.

Generic Cauchy, circulant, Hankel, and Toeplitz constructors use explicit exact rational parameter vectors accepted by the source APIs.
Random graph families and graph transformations use deterministic seed `20260802`; the original input dimension and seed are stored in
`original_id`. The graph Laplacian uses the exact unnormalized source option because normalized Laplacians can require square roots.

Special exact constructions were checked as follows:

- `cross` uses rational 3-4-5 blocks and satisfies `X^2 = I` exactly.
- `block_householder` uses the documented integer `Z` and an exact rational Gram inverse and satisfies `P^2 = I` exactly.
- `symmetric_stochastic` uses the exact Soules outer-product formula. Values for dimensions 1-15 were independently compared with the
  Anymatrix floating Givens implementation; the largest discrepancy was approximately `2.22e-16`.
- `strakos` uses exact parameters `p=1/2`, `lambda_min=1/10`, and `lambda_max=100`.
- `kms` uses exact `rho=1/2`.
- `rosser` covers supported powers of two with integer parameter pairs `(2,1)`, `(3,1)`, and `(3,2)`.
- `poisson` is included only at square output dimensions no greater than 63.

## Archived Source Matrices

All 594 in-range Anymatrix Hadamard text files were inspected. Sixteen are exactly symmetric and satisfy `H H^T = nI` using exact
integer arithmetic. Their dimensions are 1, 2, 4, 8, 16, 16, 20, 20, 20, 20, 28, 32, 36, 40, 56, and 60. Sampling retains the
required property-band and family representatives rather than all 16.

Within Anymatrix NESSIE, Benguela is an exact symmetric binary matrix of dimension 29 and is eligible. Hexgrid is an exact scalable
symmetric adjacency family. The car/train correlation matrix is excluded because its intended entries come from floating statistical
calculations rather than an exact rational source formula; `spl0708b` is asymmetric; other NESSIE matrices are rectangular,
asymmetric, or exceed dimension 63.

## Exclusions

The following are excluded rather than approximated:

- nonsymmetric or rectangular generators;
- dimensions above 63;
- irrational or transcendental formulas, including square-root, trigonometric, exponential, logarithmic, or pi-dependent entries;
- random-real generators whose output is not exact rational source data;
- normalized graph Laplacians involving square roots;
- families that become symmetric only through forced symmetrization.

Notable excluded families include the symmetric Clement mode because it uses square roots; Randcorr, Wathen, and Oscillate because they
use random-real orthogonal or factorization paths; Prolate because it uses sine and pi; and regularization families Foxgood, Gravity,
Blur, Shaw, and Ursell because their formulas involve square roots, fractional powers, exponentials, trigonometric constants, or
logarithms.

## Result

| Measure | Result |
| --- | ---: |
| Distinct eligible pool | 2,678 |
| Eligible families | 51 |
| Property classifications | 34 |
| Populated property-band strata | 166 |
| Selected distinct matrices | 484 |
| Exact matches already in database | 6 |
| Source rows proposed for insertion | 478 |
| Additional exact duplicate after circular normalization | 1 |
| Additional positive-scale equivalents of older matrices | 5 |
| Additional dimension-one equivalents of older matrix 314 | 2 |
| Diagonal Strakos rows omitted from the benchmark corpus | 20 |
| Retained new matrices | 450 |
| Retained generator families | 50 |
| New compact circular matrices | 78 |

The six reused rows are IDs 48, 49, 314, 1995, 2001, and 2155. Their tags were augmented with catalogue provenance; no duplicate row
was created. Source rows occupy IDs 2210-2687 except deleted IDs 2210-2215, 2220, 2254, and the 20 Strakos rows 2244, 2255, 2264,
2287, 2317, 2327, 2336, 2356, 2390, 2419, 2449, 2458, 2484, 2552, 2580, 2606, 2620, 2630, 2652, and 2659. ID 2215 became an exact
duplicate after circular normalization. IDs 2212-2214, 2220, and 2254 are positive rational multiples of older rows; IDs 2210-2211
were removed by the dimension-one rule in favor of ID 314. The Strakos family was later omitted because all 20 selected rows are
diagonal and the benchmark corpus retains only three diagonal representatives. Fifty-three full-storage circulant rows were
normalized by subtracting their exact common diagonal value from every entry; each description records that value. No
generator-catalogue dimension-one or Strakos row remains.

| Dimension band | New rows | New compact circular rows |
| --- | ---: | ---: |
| 1-8 | 69 | 14 |
| 9-16 | 90 | 18 |
| 17-25 | 86 | 15 |
| 26-44 | 104 | 17 |
| 45-63 | 101 | 14 |

In the original selected source set, two populated strata cannot supply three distinct eligible matrices: `fixed_size x 45-63`
has two and `unimodular x 1-8` has one. The other four `unimodular` bands have no eligible matrices. Every other populated
property-band stratum originally contained three selected representatives.

## Validation

After insertion, an independent database pass verified:

1. `PRAGMA integrity_check` returns `ok`.
2. `PRAGMA foreign_key_check` returns no rows.
3. Every stored token count matches either full upper-triangle or compact circular format.
4. Every compact row expands to the intended exact symmetric matrix.
5. Every one of the 1,064 database matrices is exactly symmetric.
6. All 1,064 stored value vectors are pairwise distinct up to a positive nonzero rational multiplier within each `(dimension,
   is_cs)` class. Twenty-seven groups related only by a negative multiplier remain intentionally distinct.
7. Every dimension is between 1 and 63.
8. Every row of dimension 26-63 has `size_class = "super_large"`, and no `super_large` tag remains.
9. The unique `(dimension, matrix)` index rejects duplicate stored matrices while permitting the same compact text at different
   dimensions.
10. One-second safe calibration supplies complete baselines for 305 retained new matrices, bringing the candidate table to 65,800
    rows; the timing table contains 710 rows.
11. One-second fast calibration records 310 measurements and 140 timeout sentinels for the 450 retained new rows, leaving no null
    fast calibration.
12. One-second safe calibration records 305 measurements and complete exact baselines plus 145 timeout sentinels, leaving no null
    safe calibration and no partial baseline.
13. A second circular-normalization audit finds zero conversions and zero duplicates.
14. The only retained diagonal matrices are IDs 314, 2180, and 2203; ID 314 is the sole dimension-one row and ID 2203 the sole
    compact all-zero game.

The final database contains 1,064 matrices: 713 complete baselines and 351 catalog-only rows.
