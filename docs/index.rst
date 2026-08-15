.. meta::
   :google-site-verification: y-jxZ5iHl68HnoLBfi3KcwCVvWzJZrMyrSujsfR2L5c

FracESSA
########

**Fractional ESS Analyzer**

FracESSA finds evolutionarily stable strategies (ESS) in symmetric games with rational payoffs. Equivalently, it finds the strict
local maximizers of a rational quadratic form over the probability simplex. It provides a C++ command-line program and a Python API.

.. image:: ../logo.png
   :alt: FracESSA logo showing a simplex in three-dimensional coordinates
   :width: 620px
   :align: center

.. container:: project-actions

   * `View FracESSA on GitHub <https://github.com/reinhardullrich/fracessa>`_
   * `Install pyfracessa from PyPI <https://pypi.org/project/pyfracessa/>`_

Mathematical problem
********************

FracESSA analyzes the standard quadratic optimization problem

.. math::

   \max_{x\in\Delta^n} x^\mathsf{T} A x,

where :math:`A` is a symmetric matrix with rational entries and :math:`\Delta^n` is the probability simplex. In the corresponding symmetric
partnership game, the ESS are exactly the strict local maximizers of this quadratic form.

Numerical contract
******************

Matrix entries, returned probability vectors, payoffs, candidate decisions in ``safe`` mode, and every final stability decision are
exact. FracESSA offers two explicitly selected search methods:

.. list-table::
   :header-rows: 1
   :widths: 15 35 50

   * - Method
     - Intended use
     - Guarantee
   * - ``safe``
     - Complete correctness-first analysis
     - Uses exact arithmetic for every candidate and stability decision.
   * - ``fast``
     - Faster exploratory analysis
     - Uses floating-point candidate filtering followed by exact checks. It can miss candidates and is not a completeness
       certificate.

The nonempty support space can contain ``2^n - 1`` supports. FracESSA generates supports one at a time, prunes strict supersets of
exact equilibria, and exploits circular symmetry when present. Candidate output is therefore an intermediate record of the ESS
search, not a census of every algebraic Nash equilibrium. One matrix uses one CPU core; Python can analyze independent matrices in
parallel.

Quick start
***********

Install the Python package:

.. code-block:: console

   python -m pip install pyfracessa

Run one exact analysis:

.. code-block:: python

   from pyfracessa import Matrix, RunConfig, run

   result = run(
       "safe",
       Matrix(1, "3#-1,0,0,-2,0,-3"),
       config=RunConfig(include_candidates=True),
   )

   print(result["ess_count"])
   print(result["candidates"])

Prebuilt CLI binaries are available from the `GitHub releases page <https://github.com/reinhardullrich/fracessa/releases>`_. Continue
with :doc:`getting-started` for matrix encoding, command-line use, result interpretation, and multiprocessing.

Documentation
*************

.. toctree::
   :maxdepth: 1
   :caption: Documentation

   getting-started
   algorithm
   python-api
   cpp-api

Project links
*************

* `Source code and test data <https://github.com/reinhardullrich/fracessa>`_
* `Standalone releases <https://github.com/reinhardullrich/fracessa/releases>`_
* `Python package <https://pypi.org/project/pyfracessa/>`_
* `Issue tracker <https://github.com/reinhardullrich/fracessa/issues>`_
* `GPL-3.0-or-later license <https://github.com/reinhardullrich/fracessa/blob/main/LICENSE>`_
