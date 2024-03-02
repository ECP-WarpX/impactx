.. _examples-from-array:

Initialize a beam from arrays
=============================

This example demonstrates how a beam can be initalized in ImpactX from array-like structures.
This allows various applications of interest, 
such as using a beam from a different simulation, 
initializing a beam from file,
or creating a custom distribution.
This example includes a set of utilities for transforming the beam to the fixed-s coordinates of ImpactX.



We use a 2 GeV electron beam.

The second moments of x, y, and t should be unchanged, but there is large emittance growth in the x and y phase planes.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can **only** be run with **Python**:

* **Python** script: ``python3 run_from_array.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: run_from_array.py
    :language: python3
    :caption: You can copy this file from ``examples/initialize_from_array/run_from_array.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_from_array.py``

   .. literalinclude:: analysis_from_array.py
      :language: python3
      :caption: You can copy this file from ``examples/intialize_from_array/analysis_from_array.py``.
