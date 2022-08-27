.. _examples-iotalens:

A nonlinear focusing channel based on the IOTA nonlinear lens
=============================================================

A constant focusing channel with nonlinear focusing, using a string of thin
IOTA nonlinear lens elements alternating with constant focusing elements.

We use a 2.5 MeV proton beam, corresponding to the nominal IOTA proton energy.

The two functions H (Hamiltonian) and I (the second invariant) should remain unchanged for all particles.

In this test, the initial and final values of :math:`\mu_H`, :math:`\sigma_H`, :math:`\mu_I`, :math:`\sigma_I` must agree with nominal values.


Run
---

This example can be run as a Python script (``python3 run_iotalens.py``) or with an app with an input file (``impactx input_iotalens.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_iotalens.py
          :language: python3
          :caption: You can copy this file from ``examples/iota_lens/run_iotalens.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_iotalens.in
          :language: ini
          :caption: You can copy this file from ``examples/iota_lens/input_iotalens.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_iotalens.py``

   .. literalinclude:: analysis_iotalens.py
      :language: python3
      :caption: You can copy this file from ``examples/iota_lens/analysis_iotalens.py``.
