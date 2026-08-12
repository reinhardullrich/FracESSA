Python API
##########

The public ``pyfracessa`` package calls the native C++ analyzer in-process. Inputs are small dataclasses and every result is an
ordinary Python dictionary; no conversion from a result class is required.

Execution contract
******************

Every execution function requires ``"safe"`` or ``"fast"``. ``safe`` is the complete exact route.
``fast`` can miss candidates during its binary64 prefilter, although every surviving candidate and final stability decision is
exact. See :doc:`getting-started` for method selection and matrix text formats.

``run()`` and ``run_multiprocessing()`` have the same return modes:

.. list-table::
   :header-rows: 1
   :widths: 32 28 40

   * - Input
     - Sink
     - Return value
   * - One :class:`~pyfracessa.Matrix`
     - None
     - One result dictionary.
   * - Iterable of matrices
     - None
     - A lazy iterator of result dictionaries.
   * - One matrix or an iterable
     - Provided
     - The number of result dictionaries written eagerly.

Result dictionary
*****************

Every result contains the following fields:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Field
     - Meaning
   * - ``run_id``
     - Identifier supplied by the caller or generated for the complete wrapper call.
   * - ``matrix_id``
     - Signed 64-bit identifier copied from the input matrix.
   * - ``status``
     - Integer :class:`~pyfracessa.StatusCode`: ``0`` success, ``1`` parse error, ``4`` execution error, or ``255`` internal error.
   * - ``candidate_count``
     - Number of exact equilibrium candidates found by the selected method, including circular multipliers. The count is complete
       with ``safe``; ``fast`` can miss candidates.
   * - ``ess_count``
     - Number of ESS found by the selected method, including circular multipliers. The count is complete with ``safe``.
   * - ``candidate_structure``
     - Dictionary mapping support size to candidate count at that size.
   * - ``ess_structure``
     - Dictionary mapping support size to ESS count at that size.
   * - ``elapsed_ns``
     - Native analyzer duration in nanoseconds, measured with a monotonic clock; Python scheduling and file output are excluded.
   * - ``safe_fallback``
     - Whole-matrix ``fast`` fallback reason, or ``None``. A local exact retry for one support does not set this field.
   * - ``error_message``
     - Empty on success; otherwise the parser, execution, or worker diagnostic.
   * - ``candidates``
     - Stored representative rows when ``RunConfig.include_candidates`` is true; otherwise an empty list.
   * - ``metadata``
     - The caller's original matrix metadata dictionary, or ``None``.

The possible non-null whole-matrix fallback reasons are ``"precision_span"``, ``"equilibration_invalid"``, and
``"equilibration_non_convergence"``.

Candidate dictionary
********************

Each stored candidate contains:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Field
     - Meaning
   * - ``candidate_id``
     - One-based identifier of the stored representative in deterministic search order.
   * - ``vector``
     - Exact full-dimensional probability vector as comma-separated rational text; entries outside the support are zero.
   * - ``support`` / ``support_size``
     - Arbitrary-precision integer bit mask for ``I(x)`` and its number of set bits.
   * - ``extended_support`` / ``extended_support_size``
     - Bit mask for ``J(x)``, all pure best replies to the candidate, and its number of set bits.
   * - ``multiplier``
     - Number of distinct circular rotations/reflections represented by the row, or ``None`` for an ordinary candidate.
   * - ``is_ess``
     - Whether exact stability checking accepted the candidate as an ESS.
   * - ``stability``
     - Machine-readable reason code for the stability decision.
   * - ``payoff``
     - Exact rational payoff as text.
   * - ``payoff_dbl``
     - Lossy binary64 approximation provided only for convenience.

Counts and structures include every candidate found by the selected method even when representative rows are disabled. For large
searches, leave candidates disabled unless their individual vectors are needed; candidate payloads can dominate memory and output
volume.

Core types
**********

.. autoclass:: pyfracessa.StatusCode
   :members:

.. autoclass:: pyfracessa.Matrix
   :members:

.. autoclass:: pyfracessa.RunConfig
   :members:

.. autoclass:: pyfracessa.MPConfig
   :members:

Analyzer
********

.. autofunction:: pyfracessa.new_run_id

.. autofunction:: pyfracessa.compute_matrix

.. autofunction:: pyfracessa.run

.. autofunction:: pyfracessa.run_multiprocessing

Multiprocessing notes
*********************

One worker analyzes one matrix at a time; parallelism is across matrices, and results are yielded in completion order. Submission is
bounded by ``min(queue_maxsize, workers * prefetch_per_worker)``. Closing an unfinished iterator terminates its workers. Scripts using
the default ``spawn`` method must protect their entry point with ``if __name__ == "__main__":``. Native logging is sequential-only;
enabling it for a multiprocessing call raises ``ValueError`` before workers are started.

Input and output
****************

.. autofunction:: pyfracessa.load_matrices_from_json

.. autofunction:: pyfracessa.create_sink

.. autoclass:: pyfracessa.CsvSink
   :members:

.. autoclass:: pyfracessa.JsonSink
   :members:

.. autoclass:: pyfracessa.ParquetSink
   :members:

All sinks create separate summary, candidate, and metadata outputs without overwriting existing paths. A failed construction, run,
write, or finalization removes only files created by that attempt and re-raises the original error. CSV and JSON preserve
arbitrary-width support integers. The current Parquet candidate schema is ``uint64`` and rejects supports above ``2^64 - 1`` rather
than truncating them.
