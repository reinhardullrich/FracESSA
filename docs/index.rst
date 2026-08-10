FracESSA Documentation
######################

FracESSA finds evolutionarily stable strategies (ESS) in symmetric games with rational payoffs. Equivalently, it finds the strict
local maximizers of a quadratic form over the probability simplex. Matrix entries, returned strategy vectors, payoffs, candidate
decisions in ``safe`` mode, and every final stability decision are exact.

Use this manual in the following order:

* :doc:`getting-started` explains installation, matrix input, method selection, results, command-line use, and Python batches.
* :doc:`algorithm` explains the mathematical problem and search pipeline from support generation through exact stability.
* :doc:`python-api` documents the Python types, execution functions, result dictionaries, multiprocessing, and file sinks.
* :doc:`cpp-api` documents the native analyzer, candidate records, and validating matrix parser.

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   getting-started
   algorithm
   python-api
   cpp-api

The nonempty support space can contain ``2^n - 1`` supports, so representable dimensions are not necessarily practical. FracESSA
generates supports one at a time, prunes strict supersets of exact equilibria, and exploits circular symmetry when present. One
matrix uses one CPU core; the Python interface can analyze independent matrices in parallel.

The `project homepage <https://reinhardullrich.github.io/fracessa/>`_ provides downloads and a short introduction. The
`GitHub repository <https://github.com/reinhardullrich/fracessa>`_ contains the complete source and test data.
