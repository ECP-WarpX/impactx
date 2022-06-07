.. _examples-triplet:

Optimized Triplet
=================

Optimization of focusing parameters for a quadrupole triplet.
A 2 GeV electron beam is strongly focused from lattice initial parameters:

:math:`\beta_x = \beta_y = 40` m
:math:`\alpha_x = -\alpha_y = 2`

to final lattice parameters:

:math:`\beta_x = \beta_y = 0.55` m
:math:`\alpha_x = \alpha_y = 0`

resulting in a reduction by a factor of 8.5 in the horizontal and vertical beam sizes.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_triplet.py`` or
* ImpactX **executable** using an input file: ``impactx input_triplet.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_triplet.py
          :language: python3
          :caption: You can copy this file from ``examples/optimize_triplet/run_triplet.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_triplet.in
          :language: ini
          :caption: You can copy this file from ``examples/optimize_triplet/input_triplet.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_triplet.py``

   .. literalinclude:: analysis_triplet.py
      :language: python3
      :caption: You can copy this file from ``examples/optimize_triplet/analysis_triplet.py``.
