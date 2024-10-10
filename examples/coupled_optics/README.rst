.. _examples-coupled-optics:

Coupled Optics
==============

This is a lattice illustrating fully coupled 6D transport.  It is obtained from the example "dogleg" by adding a solenoid after the first bending dipole.
The solenoid is identical to that found in the example "solenoid".

Its primary purpose is to benchmark the calculation of the three beam eigenemittances (mode emittances).

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must
agree with nominal values.

In addition, the initial and final values of :math:`emittance_1`, :math:`emittance_2`, :math:`emittance_3` must coincide.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_coupled_optics.py`` or
* ImpactX **executable** using an input file: ``impactx input_coupled_optics.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_coupled_optics.py
          :language: python3
          :caption: You can copy this file from ``examples/coupled_optics/run_coupled_optics.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_coupled_optics.in
          :language: ini
          :caption: You can copy this file from ``examples/coupled_optics/input_coupled_optics.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_coupled_optics.py``

   .. literalinclude:: analysis_coupled_optics.py
      :language: python3
      :caption: You can copy this file from ``examples/coupled_optics/analysis_coupled_optics.py``.
