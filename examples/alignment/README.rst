.. _examples-alignment:

Quadrupole with Alignment Errors
================================

A 2 GeV proton beam propagates through a single quadrupole with 3 mm horizontal misalignment and 30 degree rotation error.

The first and second moments of the particle distribution before and after the quadrupole should coincide with
analytical predictions, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\mu_x`, :math:`\mu_y`, :math:`\sigma_x`, :math:`\sigma_y`,
:math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_alignment.py`` or
* ImpactX **executable** using an input file: ``impactx input_alignment.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_alignment.py
          :language: python3
          :caption: You can copy this file from ``examples/alignment/run_alignment.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_alignment.in
          :language: ini
          :caption: You can copy this file from ``examples/alignment/input_alignment.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_alignment.py``

   .. literalinclude:: analysis_alignment.py
      :language: python3
      :caption: You can copy this file from ``examples/alignment/analysis_alignment.py``.
