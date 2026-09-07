.. _examples-GSI-SIS18-lattice:

The "bare" linear lattice of the Fermilab IOTA storage ring
===========================================================

The lattice of the SIS18 storage ring at GSI, in the configuration used for the following benchmark:

https://web-docs.gsi.de/~giuliano/research_activity/trapping_benchmarking/main.html

The second moments of the particle distribution after a single turn should coincide with the initial section moments of the particle distribution, to within the level expected due to numerical particle noise.
The example runs 5 turns.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_gsilattice.py`` or
* ImpactX **executable** using an input file: ``impactx input_gsilattice.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_gsilattice.py
          :language: python3
          :caption: You can copy this file from ``examples/gsi_lattice/run_gsilattice.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_gsilattice.in
          :language: ini
          :caption: You can copy this file from ``examples/gsi_lattice/input_gsilattice.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_gsilattice.py``

   .. literalinclude:: analysis_gsilattice.py
      :language: python3
      :caption: You can copy this file from ``examples/gsi_lattice/analysis_gsilattice.py``.
