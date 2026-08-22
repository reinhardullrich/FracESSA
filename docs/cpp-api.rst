C++ API
#######

The public native interface consists of the synchronous analyzer and its candidate records. The CLI and Python binding parse input
with the embedded coposit parser and pass its exact integer-scaled result directly to the analyzer. Direct C++ callers receive
exceptions from invalid method names, malformed input, allocation failures, or arithmetic failures; the CLI and Python binding
translate them into their public status formats.

Core types
**********

.. doxygenenum:: fracessa::search_method
   :project: FracESSA

``fracessa::candidate`` stores one exact symmetric Nash equilibrium. Its support handles use the matrix-wide support context, so the
same type represents every supported dimension. ``is_ess`` and ``stability`` record the later exact stability decision.

.. doxygenclass:: fracessa::candidate
   :members:
   :project: FracESSA

Analyzer
********

Constructing ``fracessa::analyzer`` runs the configured analysis synchronously. The analyzer creates one matrix-wide support context
and uses the same public search path for every dimension.

The public count and structure fields are always populated and include circular multipliers. They contain every candidate found by
the ESS search; only ``safe`` guarantees a complete ESS result. They are not an enumeration of every degenerate or strict-superset
Nash equilibrium. The ``candidates_`` vector is populated only when ``with_candidates`` was true, and omitting rows does not change
the counts. A non-``none`` ``safe_fallback_`` records only a whole-matrix switch from ``fast`` to exact candidate search.

.. doxygenfunction:: fracessa::parse_search_method
   :project: FracESSA

.. doxygenclass:: fracessa::analyzer
   :members:
   :project: FracESSA

Input and output
****************

coposit's ``coposit::parsers::matrix_parser`` accepts exact integers, fractions, finite decimals, and scientific notation in either
the full upper-triangular layout or the compact circular layout. The latter begins with the common diagonal and continues with one
value per positive circular distance. A fraction has an optional sign before its numerator
and a positive integer denominator, so ``-1/2`` is valid and ``1/-2`` is not. The parser also accepts symmetric Matrix Market input
and returns ``coposit::parsers::parsed_matrix``: one integer matrix, its common positive denominator, and the compact-circular marker
consumed by ``fracessa::analyzer``.

``matrix_parser::parse()`` reads matrix **contents**. Use ``matrix_parser::parse_file()`` when a direct C++ caller has a path. The
compact circular layout is explained in :doc:`getting-started`.

Minimal native example
**********************

The analyzer takes ownership of the parsed matrix:

.. code-block:: cpp

   #include <coposit/parsers/matrix_parser.hpp>
   #include <fracessa/fracessa.hpp>

   #include <iostream>
   #include <utility>

   int main()
   {
       auto game = coposit::parsers::matrix_parser::parse("3#-1,0,0,-2,0,-3");
       fracessa::analyzer result(fracessa::search_method::safe, std::move(game), true);
       std::cout << result.ess_count_ << '\n';
   }

Construction blocks until the complete analysis finishes. Pass ``false`` instead of ``true`` for the third argument when only counts
and support-size structures are needed; this avoids retaining individual candidate vectors without changing the computation.
