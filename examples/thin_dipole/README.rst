.. _examples-thin-dipole:

Thin Dipole
===========

This test involves tracking a 5 MeV proton beam through a 90 degree sector bend, using two different methods:

#. Using a sequence of drifts and thin kicks as described by:
G. Ripken and F. Schmidt, "A Symplectic Six-Dimensional Thin-Lens Formalism for Tracking," CERN/SL/95-12 (1995).

#. Using the exact nonlinear transfer map for a sector bend as described by:
D. Bruhwiler et al, "Symplectic Propagation of the Map, Tangent Map and Tangent Map Derivative through
Quadrupole and Combined-Function Dipole Magnets without Truncation," THP41C, EPAC98, pp. 1171-1173 (1998).

To compare the two models, the beam is tracked using 1), followed by the inverse of 2).
The beam moments and emittances should coincide with those of the initial distribution.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_thin_dipole.py`` or
* ImpactX **executable** using an input file: ``impactx input_thin_dipole.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_thin_dipole.py
          :language: python3
          :caption: You can copy this file from ``examples/thin_dipole/run_thin_dipole.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_thin_dipole.in
          :language: ini
          :caption: You can copy this file from ``examples/thin_dipole/input_thin_dipole.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_thin_dipole.py``

   .. literalinclude:: analysis_thin_dipole.py
      :language: python3
      :caption: You can copy this file from ``examples/thin_dipole/analysis_thin_dipole.py``.


Visualize
---------

.. note::

   TODO :)
