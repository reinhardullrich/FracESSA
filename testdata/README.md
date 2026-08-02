# FracESSA Test Data

`fracessa_testdata.sqlite3` is the canonical store for exact test matrices,
complete expected candidate results where available, and timing data.

The current snapshot contains 627 pairwise-distinct exact matrices. The 98
analyzed matrices have 49,388 stored candidate representatives whose
multipliers represent 86,383 candidates and 83,462 ESS; the other 529 catalog
rows deliberately have null baseline fields because exhaustive analysis has
not been run.

It contains each distinct matrix from Tables 1 and 2 of the
Bomze-Schachinger-Ullrich ESS-growth paper exactly once. IDs 18 and 26 hold the
exact published Table 1 matrices that replaced same-property alternatives;
IDs 80-90 hold the previously missing Table 2 base and constructed matrices.
Redundant alternatives formerly at IDs 12 and 21 were removed. IDs 56-66 are
staged complete-multipartite many-ESS benchmark matrices. IDs 67-79 are
deterministic random-integer coverage matrices; together with the existing rows,
every dimension from 2 through 25 has at least one circular and one non-circular
matrix. IDs 45-47 preserve the unsafe-filter, LU-boundary, and failed-proof
verified-search regressions. IDs 91-117 are SuiteSparse Matrix Collection
imports, IDs 118-269 are QAPLIB imports, and IDs 270-280 are TSPLIB95 imports.
IDs 281-313 are Biq Mac Library imports, and IDs 314-319 are Magma Hadamard
database imports. Thirty QPLIB objective imports retain their original IDs
between 320 and 1925; the gaps are the deliberately removed constraint-only
rows. IDs 1926-1982 are OR-Library imports, IDs 1983-1990 are KONECT imports,
IDs 1991-2137 are the House of Graphs stratified sample, IDs 2138-2176 are the
Network Data Repository stratified sample, and IDs 2177-2206 are SDPLIB `F0`
objective-matrix imports. No complete SQLite matrix suite is currently wired
into CTest.

The timing table contains one CPU-2 persistent-Pybind median session with a
one-second target and 507 rows over the 87 matrices that preceded the
SuiteSparse import. Catalog-only rows are excluded automatically. Current
unsafe, verified, and exact use one Release/native/LTO binary at algorithm
revision `34e003168607`; historical default, very unsafe and the retained
Werner modes use separate historical adapters. Historical very unsafe
mismatches IDs 38-39, current unsafe mismatches IDs 45-47, and current verified
and exact match all 87 matrices. Werner exact has 72 completed rows.

## SuiteSparse Imports

The SuiteSparse Matrix Collection import uses the following reproducible rule:

- The matrix is square, real, exactly symmetric after Matrix Market expansion,
  and has dimension at most 63.
- Matrix Market `pattern` entries become 1. Integer entries remain integers.
  Every finite decimal or scientific-notation `real` token is converted to its
  exact reduced fraction; this preserves the downloaded file exactly but does
  not claim that an underlying physical model was symbolically rational.
- The dense exact matrix is stored in FracESSA upper-triangle format, or compact
  circular format only when exact circular symmetry and a zero diagonal are
  present.
- Exact dense rational duplicates are not imported. SuiteSparse ID 2758,
  `Mycielski/mycielskian2`, was skipped because it equals existing matrix ID 1.
- Dimensions 26-63 are retained but tagged `"super_large"`.

The official 2,904-entry index contains 56 real square matrices through order
63. Exact parsing accepts 28 symmetric files and rejects the other 28; ID 2758
duplicates existing matrix ID 1, leaving 27 imports at IDs 91-117. IDs 91-101
retain complete exact candidate baselines from the earlier order-30 import,
while IDs 102-117 are catalog-only. Source page, collection ID, original Matrix
Market field type, group, title, date, and author are retained in metadata.

## NIST Matrix Market Audit

The finite NIST Matrix Market search over every symmetric category and
`0 < n < 64` returns ten names. Nine Matrix Market files remain downloadable;
all are exactly symmetric. Eight are exact SuiteSparse duplicates. `BFW62B`
has the same source and sparsity pattern as SuiteSparse `Bai/bfwb62`, with only
eight last-digit differences of at most `1e-21`, so it is retained as rounded
alternate provenance. NIST withdrew the incorrect `LAP 25` Matrix Market
assembly. The earlier `FIDAP005` audit likewise identifies a lower-precision
version of SuiteSparse `FIDAP/ex5` rather than a distinct game.

NIST therefore adds no matrix row. Corresponding SuiteSparse rows carry
`"nist_matrix_market"` and retain the NIST page in `description`. Parameterized
generators remain outside this finite downloadable-file audit.

## QAPLIB Imports

The official QAPLIB `qapdata` archive contains 136 fixed quadratic-assignment
instances. Each file stores one dimension followed by two dense integer
matrices, A and B. The import tests A and B independently, retains only exact
symmetric matrices with dimensions 1-63, and globally deduplicates their dense
values before choosing compact circular or upper-triangle FracESSA storage.

Of 109 in-range instance files, 182 matrix occurrences are symmetric. Thirty
are repeated matrices within QAPLIB, leaving 152 exact distinct imports at IDs
118-269; one is circular-symmetric. These rows are catalog-only. The source is
the QAPLIB dataset by Burkard, Cela, Karisch, Rendl, Anjos, and Hahn,
[DOI 10.7488/ds/3428](https://doi.org/10.7488/ds/3428), licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The database stores
the primary `instance:role` identifier and lists any identical alternate source
occurrences in `description`.

## TSPLIB95 Imports

The official [TSPLIB95 symmetric-TSP archive](https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/ALL_tsp.tar.gz)
contains 111 problem files. Seventeen declare `EDGE_WEIGHT_TYPE: EXPLICIT`; six
have dimensions above FracESSA's limit, leaving 11 exact symmetric integer
matrices with dimensions 17-58. Coordinate-derived distances, tours, and the
separate asymmetric problem categories are not imported.

The retained `FULL_MATRIX`, `LOWER_DIAG_ROW`, and `UPPER_ROW` representations
were expanded exactly and checked edge-for-edge against TSPLIB95's independent
official XML files. All 11 are non-circular and distinct from each other and
from every existing dense exact database matrix. They are catalog-only rows at
IDs 270-280; dimensions 26-63 carry the `"super_large"` tag.

## Biq Mac Library Imports

The official [Biq Mac Library archive](https://biqmac.aau.at/library/tar_files/biqmac_all.tar.gz)
contains 468 files representing 343 logical binary-quadratic and Max-Cut
instances. Its 125 dense/sparse pairs match exactly. Applying the dimension-63
limit leaves 33 matrices: 10 Beasley binary-quadratic instances, 13
Glover-Kochenberger-Alidaee binary-quadratic instances, and 10 Rudy Max-Cut
graphs. The other 310 logical instances are too large.

All retained entries are exact integers. Sparse binary-quadratic entries are
expanded under the library's symmetric-Q contract; each Max-Cut graph becomes
its symmetric zero-diagonal weighted adjacency matrix. The 46 individual files
behind the retained rows match their aggregate-archive copies byte-for-byte.
All 33 matrices are non-circular and globally distinct. They are catalog-only
rows at IDs 281-313; dimensions 26-63 carry the `"super_large"` tag. The audited
archive SHA-256 is `887ed2a8187fff2cf941d3c6aad3953ffd9904700dd86710fc2bd09736670e5a`.

## Magma Hadamard Imports

The official Magma downloads provide separate archives for
[ordinary Hadamard matrices](https://magma.maths.usyd.edu.au/magma/download/db/hadamard.tar.gz)
and [skew-Hadamard matrices](https://magma.maths.usyd.edu.au/magma/download/db/hadamard_skew.tar.gz).
Both downloads match Magma's published MD5 and SHA-1 values. Their compact
binary representation was decoded independently and checked against the exact
degree-16 representative printed in the Magma handbook.

The ordinary database contains 5,391 inequivalent representatives, of which
4,474 have degree at most 63. Every in-range matrix satisfies `H H^T = nI`, but
only six are themselves symmetric: representative 1 at each of degrees 1, 2,
4, 8, 16, and 32. All 638 matrices in the separate skew database have degrees
36, 44, or 52; all satisfy the skew-Hadamard property and none is symmetric.
The six retained exact `+/-1` matrices are non-circular and globally distinct.
They are catalog-only rows at IDs 314-319; degree 32 carries the
`"super_large"` tag.

The audited SHA-256 values are
`69930089fe46dd59cb0c48c73e1cfd03928b2e25958b1ce22de7a9f647e0cad7`
for the ordinary archive and
`1aa7f7fef8e541c7078ed89431a42a1814a786d74fb0f0d777f06babded5f210`
for the skew archive.

## QPLIB Imports

The official [QPLIB](https://qplib.zib.de/) catalog contains 453 quadratic-
optimization instances and is licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The import uses the
35 canonical `.qplib` files with 1-63 variables rather than downloading the
775 MB multi-format bundle. QPLIB's defining publication is available at
[DOI 10.1007/s12532-018-0147-4](https://doi.org/10.1007/s12532-018-0147-4).

Only an explicitly stored quadratic objective is imported. QPLIB stores its
lower triangle, which is mirrored exactly to recover the symmetric Hessian.
Quadratic-constraint matrices, linear objective vectors, and generally
rectangular linear-constraint matrices are excluded. Finite decimal and
scientific-notation coefficients are converted to their exact reduced
fractions, preserving the source tokens without claiming a more precise
underlying physical value.

Thirty of the 35 files have an explicit quadratic objective. Their 30 objective
matrices are pairwise distinct, duplicate no earlier database row, and are not
circular. Every coefficient and role was independently cross-checked with
PyQPLIB 0.1.7. The rows retain their original noncontiguous IDs between 320 and
1925, are catalog-only, and carry the `"super_large"` tag at dimensions 26-63.

## OR-Library Imports

The audit covers every locally maintained problem family in J.E. Beasley's
official [OR-Library index](https://people.brunel.ac.uk/~mastjjb/jeb/info.html),
plus its still-hosted urban-transit page. OR-Library material is published
under the [MIT license](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/legal.html).
Only explicitly stored square matrices with exact rational entries, exact
symmetry, and dimensions 1-63 are eligible. Rectangular problem tables,
coordinate-derived distances, shortest-path matrices derived from edge lists,
and external collections merely linked by OR-Library are outside this import.

The retained data comprise 23 binary-quadratic Q matrices, 23 capacitated
minimum-spanning-tree cost matrices, six aircraft-separation matrices, two CAB
hub-location matrices, one portfolio correlation matrix, and two urban-transit
demand matrices. The binary-quadratic files state a maximization problem and
have the opposite sign from the corresponding Biq Mac minimization copies, so
they are distinct games rather than duplicates.

The in-range source audit found 68 symmetric matrix occurrences. Ten matrices
are repeated between `capmst1` and `capmstnew.zip`, and `portreb1` repeats the
`port1` correlation matrix, leaving 57 exact distinct imports at IDs 1926-1982.
Two additional in-range `capmstnew.zip` tables, two aircraft-separation tables,
and all six in-range corporate withholding-tax tables are asymmetric. The
urban-transit demand files `td1` and `td2` were recovered from 2011 Internet
Archive snapshots of the now-missing official files; `td3` exceeds dimension
63, while the time files contain nonnumeric `-` entries for absent links and
are not rational matrices. All 57 imports are non-circular, globally new, and
catalog-only; dimensions 26-63 carry the `"super_large"` tag.

## COMPl_e_ib Audit

The official [COMPl_e_ib 1.1](https://www.compleib.de/) archive defines 168
control-system benchmarks. Its benchmark state matrix is `A`; the other output
arrays describe input, output, noise, and weighting channels and are not
independent benchmark matrices. Of the 168 state matrices, 57 exceed
FracESSA's dimension-63 limit. All 111 in-range state matrices were constructed
from the archive's `COMPleib.m` and `.mat` files and compared exactly with their
transposes. None is symmetric, so COMPl_e_ib contributes no catalog row.

Square identity, zero, and weighting arrays synthesized while assembling a
control problem are deliberately not imported as separate games.

## SLICOT Model-Reduction Audit

The official [SLICOT model-reduction collection](https://www.slicot.org/20-site/126-benchmark-examples-for-model-reduction)
contains 18 linear-system benchmarks. Seventeen have orders from 84 to 10,913
and exceed FracESSA's dimension limit. The only in-range case is the order-48
building model `build.mat`; its stored state matrix `A` is not symmetric. The
complete MAT file is byte-identical to COMPl_e_ib's already-audited `lah.mat`,
so the collection contributes no catalog row.

## KONECT Imports

The official [KONECT network collection](https://konect.cc/networks/) currently
lists 1,326 datasets. The in-range audit downloads all 23 available unipartite
networks with at most 63 vertices and reconstructs their adjacency matrices
exactly from KONECT's whitespace-separated edge files. Unweighted edges become
1; signed and positive integer weights are retained exactly; repeated directed
edges, where present, are summed at their matrix position. Seven in-range
bipartite networks are excluded because their native adjacency tables are
rectangular rather than FracESSA game matrices.

All nine undirected files are symmetric. One of the 14 directed files,
`moreno_taro`, is also exactly symmetric because every directed edge has its
reciprocal; the other 13 directed matrices fail exact symmetry. The Dolphins
and Zachary karate-club matrices exactly duplicate existing SuiteSparse IDs 114
and 115, which now retain KONECT as alternate provenance. The other eight
matrices are globally new, non-circular, catalog-only rows at IDs 1983-1990.
Dimensions 26-63 carry the `"super_large"` tag. KONECT publishes the files and
their source citations but does not state one collection-wide dataset license,
so this catalog makes no broader licensing claim.

## House of Graphs Sample

The official [House of Graphs](https://houseofgraphs.org/) contained 29,139 graphs on August 2, 2026. Its complete
order-1-through-63 query returned 23,988 graphs. Every canonical graph6 string was decoded as an exact, unweighted,
undirected, simple graph and checked against the API adjacency list before its zero-diagonal `0/1` adjacency matrix was
considered for import.

The sample crosses five dimension bands (`1-8`, `9-16`, `17-25`, `26-44`, and `45-63`) with ten categories: acyclic,
connected bipartite cyclic, connected planar non-bipartite, connected nonplanar, disconnected cyclic, regular, dense
(density at least 0.5), vertex-transitive, asymmetric, and an unrestricted control. For each populated stratum, the three
lowest SHA-256 ranks of `20260802|dimension_band|category|graph_id` were selected. This fixed rule makes the nominally
random sample reproducible without storing or depending on API result order.

All 50 strata were populated. The dimension-45-through-63 disconnected-cyclic stratum contained only one graph, and House
of Graphs ID 21044 was selected independently by two strata. Removing that overlap leaves 147 unique matrices: 29, 30, 30,
30, and 28 from the five dimension bands. None exactly duplicates an existing dense matrix; eight admit compact circular
storage. The catalog-only rows occupy IDs 1991-2137, retain their House of Graphs ID, canonical graph6 string, source name
when available, selected strata, source page, and seed, and tag dimensions 26-63 as `"super_large"`.

House of Graphs documents canonicalization, duplicate rejection, and downloadable formats for further personal use, but
the audit found no collection-wide data-license statement. The catalog therefore records provenance without making a
broader licensing claim.

## Network Data Repository Sample

The current [Network Data Repository](https://networkrepository.com/network-data.php) category indexes expose 6,628 graph
rows. Of these, 1,241 report dimensions 1-63 across 15 categories. The deterministic sample crosses each category with five
dimension bands (`1-8`, `9-16`, `17-25`, `26-44`, and `45-63`) and ranks candidates by SHA-256 with seed `20260802`. Ranked
candidates are examined until at most three globally new, directly represented symmetric matrices are found per populated
category/band.

The source archive, not merely the index metadata, determines eligibility. A Matrix Market file must be square and either
declare symmetric storage or satisfy exact symmetry after full expansion. An edge list must be explicitly undirected in its
metadata or comments, or contain exact reciprocal edges. Arbitrary external vertex IDs are replaced by their sorted dense
order. Integer, decimal, and fractional weights remain exact. Rectangular matrices, unsymmetric matrices, temporal edge
streams, bipartite tables without a square adjacency matrix, malformed files, and matrices requiring forced symmetrization
are excluded. Exact duplicates already in the 558-row pre-import database are skipped without changing the existing row.

The resulting 39 catalog-only imports occupy IDs 2138-2176: 15 animal-social, 15 cheminformatics, six protein, two DIMACS,
and one biological matrix. Their dimension-band counts are 6, 6, 7, 10, and 10 respectively. Twenty carry the
`"super_large"` tag, 11 have non-unit weights, and one admits compact circular storage. All 39 source archives round-trip to
their stored exact matrices, and none duplicates another database row. The repository's
[data policy](https://networkrepository.com/policy.php) states a Creative Commons Attribution-ShareAlike license without
naming a version; each row retains its dataset page for source-specific attribution.

## SDPLIB Imports

The [SDPLIB 1.2 mirror](https://github.com/vsdp/SDPLIB) contains 92 semidefinite-programming test problems in SDPA sparse
block-diagonal format. Thirty problems have an aggregate matrix dimension from 1 through 63. For each eligible problem,
the import retains only matrix number zero, `F0`, as the objective coefficient matrix. The 1,799 `F1...Fm` constraint
matrices are deliberately excluded, and the source blocks are expanded into one complete matrix rather than cataloged as
independent submatrices.

IDs 2177-2206 are the resulting 30 catalog-only objective matrices: four control, 15 H-infinity, two infeasible-dual, two
infeasible-primal, three quadratic-assignment, one theta, and three truss problems. Exact parsing converts finite decimal
and scientific-notation tokens to reduced fractions. The matrices are pairwise distinct and duplicate no earlier database
row; 15 carry the `"super_large"` tag. `theta1` is circulant but has a nonzero diagonal, so it uses the full upper-triangle
representation rather than FracESSA's compact zero-diagonal circular format. Each row links to its exact source file. The
current GitHub mirror declares GPL-3.0; this catalog makes no broader licensing claim.

## Tables

### `matrices`

Each row stores one exact matrix input and its summary:

- `matrix_id`: stable FracESSA verification ID.
- `dimension`: number of strategies, from 1 through 63.
- `size_class`: `small` for dimensions 1-8, `medium` for 9-16, and
  `large` for 17-63.
- `is_cs`: 1 for compact circular-symmetric input, otherwise 0.
- `matrix`: the exact comma-separated input values, without the `n#` prefix.
- `candidate_count` and `ess_count`: complete weighted baseline counts, or null
  together when the matrix is cataloged but not analyzed.
- `candidate_structure`: JSON object mapping support size to weighted candidate
  count.
- `ess_structure`: JSON object mapping support size to weighted ESS count. Both
  structure fields are null exactly when both count fields are null.
- `origin`: where the matrix came from and why it was retained.
- `tags`: JSON array of short qualitative categories such as
  `"numerical_edge"`, `"support_frontier"`, or `"super_large"`.
- `name`: short human-readable matrix name.
- `family` and `subfamily`: broad and narrow matrix classifications used for
  selecting related fixtures.
- `description`: fuller provenance and purpose text.
- `source_url`: original website or DOI when an external source is known.
- `original_format`: source representation or construction method.
- `original_id`: identifier used by the source, when one exists.
- `created_at`: first known project use or current-row creation as
  `YYYY-MM-DD`.

The stable `matrix_id` remains the program-facing identifier. `origin` preserves
the existing human provenance text, while `source_url` provides a
machine-readable external origin without overloading that text. A null
`source_url` or `original_id` means that the matrix was constructed inside the
project or that no external source was recorded; it does not mean the source was
searched exhaustively.

Legacy rows are marked by `family = "historical"`. When their exact insertion
day is unknown, `created_at` uses January 1 of the known year as an explicit
placeholder; for example, `2014-01-01` means only "legacy matrix known from
2014." Its purpose is to distinguish long-standing matrices from cases added to
the current suite in 2026, not to assert an exact historical day.

For example, `{"1":8,"4":2}` means eight support-one results and two
support-four results. Empty ESS structure is stored as `{}`.

### `candidates`

Each row mirrors one candidate output row. A circular row stores only the
smallest support in its rotation/reflection orbit, and its non-null `multiplier`
is the number of distinct supports represented. A non-circular row has a null
multiplier. Exact fractions and vectors remain text; `payoff_double` is also
text so database reads cannot alter formatting. `(matrix_id, candidate_id)` is
the primary key, and a support may occur only once for a matrix.

Fixed facts already represented by columns, including size, circular symmetry,
counts, and support-size structures, are not duplicated in `tags`.

### `timings`

Each row is one sequential analyzer timing measurement for one matrix. A
session may contain several builds, but each build is measured by a separate
runner invocation. Rows record the machine and pinned CPU, human build label,
moving source reference such as `main`, immutable revision, binary SHA-256,
Pybind or CLI backend, numerical mode, target and measured wall durations,
iteration count, median native duration in nanoseconds, observed ESS count,
and an optional comment.

The observed count remains separate from the expected count in `matrices`, so a
report can expose unsafe-mode mismatches without hiding or rejecting their
timings. The report derives the Bomze-Schachinger-Ullrich exponential-growth
lower bound `expected_ess ** (1 / dimension)` and prints it with dimension and
circularity; this value is not stored in the database. Old CLI builds simply
have no `safe` rows.

## Scope

`python -m fracessa.timing` reads matrices from this database and writes timing
observations back to `timings`. It is deliberately a sequential, Linux
CPU-affinity runner, not part of the multiprocessing wrapper. One Pybind process
stays open for all selected modes and matrices in a build. The first returned
C++ duration chooses `ceil(target / duration)` total samples and remains part
of the sample; a duration at or above the target chooses one run. The stored
result is the median returned `elapsed_ns`. Wall time is recorded as metadata
but does not choose the sample count or result. The tool also supports legacy
CLI timers whose output unit is supplied on the command line, but those
fresh-process rows must not be mixed with persistent-Pybind microbenchmarks. No
active matrix-verification runner is wired into CTest yet.

The schema is defined in `schema.sql`. The C++ runtime does not read this
database; the timing tool uses Python's standard `sqlite3` module.

## Integrity

Basic database checks are:

```bash
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA integrity_check;'
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA foreign_key_check;'
```
