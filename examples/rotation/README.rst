.. _examples-rotation:

Drift using PROT
========================

A drift that takes place in a rotated frame, using initial and final
applications of PROT.

We use a 2 GeV electron beam.

The second moments of x, y, and t should be nearly unchanged.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as a Python script (``python3 run_rotation.py``) or with an app with an input file (``impactx input_rotation.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_rotation.py
          :language: python3
          :caption: You can copy this file from ``examples/rotation/run_rotation.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_rotation.in
          :language: ini
          :caption: You can copy this file from ``examples/rotation/input_rotation.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_rotation.py``

   .. literalinclude:: analysis_rotation.py
      :language: python3
      :caption: You can copy this file from ``examples/rotation/analysis_rotation.py``.
