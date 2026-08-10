Getting Started
###############

Python
------

Install the published package:

.. code-block:: console

   python -m pip install pyfracessa

Run one exact analysis:

.. code-block:: python

   from pyfracessa import Matrix, RunConfig, run

   result = run("safe", Matrix(1, "2#0,1,0"), config=RunConfig())
   print(result["ess_count"])

Command line
------------

Download a binary from the `release page <https://github.com/reinhardullrich/fracessa/releases>`_, or build the CLI from source:

.. code-block:: console

   cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
   cmake --build cpp/build --parallel
   ./cpp/build/fracessa safe "3#-1,0,0,-2,0,-3"

Input uses ``dimension#values``. Values are either the upper triangle of a symmetric matrix or the compact circular-symmetric form.
Choose ``safe`` for a fully exact candidate search. ``fast`` uses binary64 candidate filtering with exact fallbacks and exact final
stability decisions, but remaining heuristic rejections mean it is not a correctness certificate.
