Getting Started
###############

FracESSA accepts exact symmetric rational matrices through either Python or the command line. Both interfaces call the same C++
analyzer and return the same mathematical counts.

Install
*******

Python 3.11 through 3.14 users can install the native extension and Python API from PyPI:

.. code-block:: console

   python -m pip install pyfracessa

Parquet output is optional:

.. code-block:: console

   python -m pip install "pyfracessa[parquet]"

Standalone command-line binaries for Linux, macOS, and Windows are available on the
`release page <https://github.com/reinhardullrich/fracessa/releases>`_. To build the CLI and Python extension from source instead:

* install a C++17 compiler, CMake 3.18 or newer, Python 3.11 or newer with development headers, GMP, MPFR, and FLINT;
* use Git with submodules enabled; and
* allow the first CMake configure to download ``argparse``, ``pybind11``, ``spdlog``, and GoogleTest.

.. code-block:: console

   git clone --recurse-submodules https://github.com/reinhardullrich/fracessa.git
   cd fracessa
   cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
   cmake --build cpp/build --parallel

The resulting command-line program is ``cpp/build/fracessa`` and the Python extension is ``cpp/build/fracessa_core.*``.

Choose a method
***************

Every call requires an explicit method; there is no default.

.. list-table::
   :header-rows: 1
   :widths: 12 35 53

   * - Method
     - Intended use
     - Numerical contract
   * - ``safe``
     - Complete, correctness-first analysis
     - Uses exact arithmetic for every candidate and stability decision.
   * - ``fast``
     - Faster analysis when completeness is not required
     - Uses a binary64 candidate filter. Every surviving candidate and every final stability decision is checked exactly, but an
       early floating-point rejection can miss a valid candidate.

``fast`` switches the complete matrix to exact candidate search when preparation detects a dangerous numerical case.
The result then records the reason in ``safe_fallback``. An inconclusive solve for only one support is also retried exactly, but that
local retry does not set ``safe_fallback``. A candidate returned by ``fast`` is therefore a genuine, exactly verified candidate; the
risk is an incomplete result because the floating-point filter can reject a valid support too early.

Encode a matrix
***************

Every input has the form:

.. code-block:: text

   dimension#comma-separated-values

Values can be integers, fractions, finite decimals, or scientific numbers, such as ``-3/5``, ``-0.6``, or ``-6e-1``. Every form is
parsed as an exact rational number. A fraction's denominator must be a positive integer; put an optional sign before the numerator,
writing ``-1/2`` rather than ``1/-2``. There are two accepted value layouts.

General symmetric matrix
------------------------

Write the upper triangle row by row. For

.. math::

   A=\begin{pmatrix}
   a_{11}&a_{12}&a_{13}\\
   a_{12}&a_{22}&a_{23}\\
   a_{13}&a_{23}&a_{33}
   \end{pmatrix},

the input is:

.. code-block:: text

   3#a11,a12,a13,a22,a23,a33

Exactly ``n(n+1)/2`` values are required.

Compact circular-symmetric matrix
---------------------------------

A circular-symmetric matrix is unchanged when all row and column indices are shifted cyclically by the same amount. Because the
matrix is also symmetric, an entry depends only on the circular distance

.. math::

   \operatorname{dist}(i,j)=\min\left(|i-j|,\,n-|i-j|\right).

For dimensions two and larger, compact input supplies one payoff for each nonzero circular distance:

.. code-block:: text

   n#c1,c2,...,c_floor(n/2)

Exactly ``floor(n/2)`` values are required. The diagonal is omitted because one common diagonal value can always be normalized to
zero without changing the problem. If the original diagonal value is :math:`d`, replace the complete matrix by

.. math::

   \widetilde A=A-d\mathbf 1\mathbf 1^\mathsf{T}.

For every :math:`x` in the simplex,

.. math::

   \widetilde A x=Ax-d\mathbf 1,
   \qquad
   x^\mathsf{T}\widetilde A x=x^\mathsf{T}Ax-d.

Thus every pure-strategy payoff and every objective value is shifted by the same constant. Best replies, candidates, local
maximizers, and ESS are unchanged. In practical terms, subtract :math:`d` from **every** original matrix entry before writing the
compact values; do not merely replace the diagonal by zero.

For example, a circular matrix with diagonal ``2`` and successive distance payoffs ``3`` and ``5`` is normalized to ``5#1,3``. That
compact input expands to a 5-by-5 matrix whose first row is ``0,1,3,3,1``. FracESSA recognizes the shorter value count automatically
and applies its circular-symmetry reductions. Dimension one and circular matrices that should not be normalized must use the full
upper-triangular layout.

Matrix files
------------

The CLI accepts a file path in place of inline input. The file can contain either ``dimension#values`` or a Matrix Market ``array``
or ``coordinate`` matrix declared ``symmetric``. Matrix Market integer, real, pattern, and complex fields are accepted; complex
entries must have zero imaginary part. Decimal and scientific values remain exact. Slash fractions are available only in
``dimension#values`` input.

Python accepts the same two **contents** in ``Matrix.matrix``; it does not interpret that string as a file path. Use
``Path.read_text()`` to load a file. Python also accepts a values-only string when the matrix metadata contains an integer
``dimension``:

.. code-block:: python

   Matrix(7, "-1,0,0,-2,0,-3", metadata={"dimension": 3})

Run one matrix from Python
**************************

.. code-block:: python

   from pyfracessa import Matrix, RunConfig, run

   matrix = Matrix(1, "3#-1,0,0,-2,0,-3")
   result = run(
       "safe",
       matrix,
       config=RunConfig(include_candidates=True),
       run_id="example",
   )

   print(result["status"])
   print(result["ess_count"])
   print(result["ess_structure"])
   print(result["candidates"])

The result is a plain dictionary. ``candidate_count`` and ``ess_count`` are the totals found by the selected method, and
``candidate_structure`` and ``ess_structure`` divide those totals by support size. These four fields are always returned. The ESS
total is complete with ``safe``; candidate totals follow FracESSA's ESS-search semantics and are not a census of every degenerate or
pruned-superset equilibrium. ``fast`` can miss valid supports. ``include_candidates=True`` also returns each stored candidate
representative; it does not change the counts or the analysis.

Check ``status`` before consuming a result. ``0`` means success, ``1`` means invalid matrix input, ``4`` means an execution error,
and ``255`` means an unexpected internal error. A nonzero status is accompanied by ``error_message``. ``elapsed_ns`` measures only
the native analyzer with a monotonic clock; Python scheduling and output writing are excluded. FracESSA deliberately imposes no
per-matrix timeout.

Run one matrix from the command line
************************************

The general form for inline input is:

.. code-block:: console

   ./cpp/build/fracessa [OPTIONS] METHOD "DIMENSION#VALUES"

For example:

.. code-block:: console

   ./cpp/build/fracessa safe --candidates "3#-1,0,0,-2,0,-3"

The CLI always writes one JSON summary line. ``--candidates`` appends semicolon-separated candidate rows. Useful options are:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Meaning
   * - ``--candidates``
     - Print stored candidate representatives after the summary.
   * - ``--fullsupport``
     - Check the full support first. If the selected candidate route finds it and exact stability accepts it as an ESS, stop;
       otherwise continue the normal search without checking it twice.
   * - ``--log``
     - Write a rotating diagnostic trace to ``log/fracessa.log`` relative to the current working directory.
   * - ``--matrixid ID``
     - Attach a signed 64-bit matrix identifier to the summary and log.

The decimal ``support`` masks in candidate rows use the least significant bit for the first matrix strategy. For example, support
``7`` has its first three bits set. On failure the CLI still writes a JSON summary to standard output, writes a diagnostic to
standard error, and exits unsuccessfully.

Run several matrices
********************

``run()`` processes matrices sequentially. ``run_multiprocessing()`` uses several worker processes, each analyzing one matrix at a
time. Results arrive in completion order, so use ``matrix_id`` rather than list position to identify them.

.. code-block:: python

   from pyfracessa import MPConfig, Matrix, run_multiprocessing

   matrices = [
       Matrix(1, "3#-1,0,0,-2,0,-3"),
       Matrix(2, "3#4,13/2,1/2,5,11/2,3"),
   ]

   if __name__ == "__main__":
       for result in run_multiprocessing("safe", matrices, mp_config=MPConfig(workers=4)):
           print(result["matrix_id"], result["ess_count"])

The ``if __name__ == "__main__":`` guard is required with the portable default start method, ``spawn``. Native logging is available
only for sequential runs. For large streams, consume the iterator continuously or pass a CSV, JSON, or Parquet sink so results do
not accumulate in memory.

Dimensions and running time
***************************

Dimensions 1 through 64 use one 64-bit support word; larger dimensions use a runtime-sized multiword mask. There is no separate
application-level dimension cap, but the dense matrix must fit in memory and the nonempty support space contains ``2^n - 1``
possibilities. Full-support searches and strongly pruned games can therefore be practical above 64 even though exhaustive searches
usually are not.

Continue with :doc:`algorithm` for the search itself, :doc:`python-api` for every Python field and return mode, or :doc:`cpp-api`
for the native interface.
