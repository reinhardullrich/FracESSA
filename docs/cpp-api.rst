C++ API
#######

The native public interface consists of the analyzer, its result rows, and the validating matrix parser. The CLI and Python binding
use the same interfaces.

Search method
-------------

.. doxygenenum:: search_method
   :project: FracESSA

.. doxygenfunction:: parse_search_method
   :project: FracESSA

Analyzer
--------

``fracessa`` is the one-word specialization of ``basic_fracessa`` used through dimension 64. ``multiword_fracessa`` is the
corresponding specialization for larger dimensions.

.. doxygenclass:: basic_fracessa
   :members:
   :project: FracESSA

Candidate result
----------------

``candidate`` and ``multiword_candidate`` are the corresponding specializations of ``basic_candidate``.

.. doxygenclass:: basic_candidate
   :members:
   :project: FracESSA

Matrix parser
-------------

.. doxygenfunction:: matrix_parser::parse_matrix_string
   :project: FracESSA
