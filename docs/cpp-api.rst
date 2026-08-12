C++ API
#######

The public native interface consists of the synchronous analyzer, its candidate records, and the validating exact matrix parser.
The CLI and Python binding use these same interfaces. Direct C++ callers receive exceptions from invalid method names, malformed
input, allocation failures, or arithmetic failures; the CLI and Python binding translate them into their public status formats.

Core types
**********

.. doxygenenum:: search_method
   :project: FracESSA

``candidate`` is the one-word specialization of ``basic_candidate`` used through dimension 64. ``multiword_candidate`` is the
corresponding specialization for larger dimensions. A candidate is already an exact symmetric Nash equilibrium; ``is_ess`` and
``stability`` record the later exact stability decision.

.. doxygenclass:: basic_candidate
   :members:
   :project: FracESSA

Analyzer
********

Constructing ``basic_fracessa`` runs the configured analysis synchronously. ``fracessa`` is its one-word specialization through
dimension 64; ``multiword_fracessa`` is the multiword specialization for larger dimensions. Select the specialization from the
parsed matrix dimension before construction, as the CLI and Pybind adapter do.

The public count and structure fields are always populated and include circular multipliers. They contain every candidate found by
the selected method; only ``safe`` guarantees a complete search. The ``candidates_`` vector is populated only when
``with_candidates`` was true, and omitting rows does not change the counts. A non-``none`` ``safe_fallback_`` records only a
whole-matrix switch from ``fast`` to exact candidate search.

.. doxygenfunction:: parse_search_method
   :project: FracESSA

.. doxygenclass:: basic_fracessa
   :members:
   :project: FracESSA

Input and output
****************

The parser accepts exact integers and integer fractions in either the full upper-triangular layout or the compact zero-diagonal
circular layout. It replaces both output arguments only after interpreting the payload and reports malformed input with
``std::invalid_argument``.

.. doxygenfunction:: matrix_parser::parse_matrix_string
   :project: FracESSA
