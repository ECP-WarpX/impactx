.. _examples-iotalattice:

The "bare" linear lattice of the Fermilab IOTA storage ring
===========================================================

The linear lattice of the IOTA storage ring, configured for operation with a 2.5 MeV proton beam.

The drift regions available for insertion of the special nonlinear magnetic element for integrable optics experiments are denoted ``dnll``.

The second moments of the particle distribution after a single turn should coincide with the initial section moments of the particle distribution, to within the level expected due to numerical particle noise.
The example runs 5 turns.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_iotalattice.py`` or
* ImpactX **executable** using an input file: ``impactx input_iotalattice.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_iotalattice.py
          :language: python3
          :caption: You can copy this file from ``examples/iota_lattice/run_iotalattice.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_iotalattice.in
          :language: ini
          :caption: You can copy this file from ``examples/iota_lattice/input_iotalattice.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_iotalattice.py``

   .. literalinclude:: analysis_iotalattice.py
      :language: python3
      :caption: You can copy this file from ``examples/iota_lattice/analysis_iotalattice.py``.


.. _examples-iotalattice-sdep:

The full nonlinear lattice of the Fermilab IOTA storage ring
=============================================================

The full nonlinear lattice of the IOTA storage ring, configured for operation with a 2.5 MeV proton beam.

The special nonlinear magnetic element for integrable optics experiments is denoted ``nll``.  To simplify analysis,
the lattice has been arranged so that nll appears as the first element in the sequence.

The two functions H (Hamiltonian) and I (the second invariant) are evaluated at the entrance to the nonlinear element.
These values should be unchanged for all particles (to within acceptable tolerance), over the specified number of periods (default 5).

In this test, the initial and final values of :math:`\mu_H`, :math:`\sigma_H`, :math:`\mu_I`, :math:`\sigma_I` must agree with nominal values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_iotalattice_sdep.py`` or
* ImpactX **executable** using an input file: ``impactx input_iotalattice_sdep.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_iotalattice_sdep.py
          :language: python3
          :caption: You can copy this file from ``examples/iota_lattice/run_iotalattice_sdep.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_iotalattice_sdep.in
          :language: ini
          :caption: You can copy this file from ``examples/iota_lattice/input_iotalattice_sdep.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_iotalattice_sdep.py``

   .. literalinclude:: analysis_iotalattice_sdep.py
      :language: python3
      :caption: You can copy this file from ``examples/iota_lattice/analysis_iotalattice_sdep.py``.
