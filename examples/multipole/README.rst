.. _examples-multipole:

Chain of thin multipoles
========================

A series of thin multipoles (quad, sext, oct) with both normal and skew coefficients.

We use a 2 GeV electron beam.

The second moments of x, y, and t should be unchanged, but there is large emittance growth in the x and y phase planes.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as a Python script (``python3 run_multipole.py``) or with an app with an input file (``impactx input_multipole.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_multipole.py
          :language: python3
          :caption: You can copy this file from ``examples/multipole/run_multipole.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_multipole.in
          :language: ini
          :caption: You can copy this file from ``examples/multipole/input_multipole.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_multipole.py``

   .. literalinclude:: analysis_multipole.py
      :language: python3
      :caption: You can copy this file from ``examples/multipole/analysis_multipole.py``.
