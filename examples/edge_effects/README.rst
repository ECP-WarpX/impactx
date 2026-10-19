.. _dipedge-nonlinear:

Error scaling test for nonlinear dipole fringe field map
========================================================

This benchmark tests the use of the nonlinear ``DipEdge`` model for integrating through a dipole fringe field.

Six distinct initial conditions are tested for a nominal proton beam with 800 MeV kinetic energy.  The values of the field integrals K0-K6 are set to the default
values, as used in:

K. Hwang and S. Y. Lee, "Dipole fringe field map for compact synchrotrons," Phys. Rev. Accel. Beams 18, 122401 (2015)

The initial conditions are chosen with increasing distance from the origin in phase space.  The value of the Lie generator is a dynamical
invariant of the ideal fringe field map.  In reality, there is an error in the final variables that scales with (x,px,y,py,t,pt,g) as degree 3.  As a result, the
Lie generator is not an exact invariant of the numerically-computed fringe field map.

In this test, the change in the Lie generator for each initial condition should coincide with its (small) nominal value.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_dipedge.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_dipedge.py
          :language: python3
          :caption: You can copy this file from ``examples/edge_effects/run_dipedge.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_dipedge.py``

   .. literalinclude:: analysis_dipedge.py
      :language: python3
      :caption: You can copy this file from ``examples/edge_effects/analysis_dipedge.py``.


.. _dipedge-linear-map:

Benchmark of linear dipole fringe field map
============================================

This benchmark tests the use of the linear ``DipEdge`` model for integrating through a dipole fringe field.

Six distinct initial conditions are tested for a nominal proton beam with 800 MeV kinetic energy.

In this test, all 36 elements of the 6x6 transport matrix must agree with predicted values obtained from the model in

K. Bane, "A First and Second-Order Matrix Theory for the Design of Beam Transport Systems and Charged Particle Spectrometers," SLAC-75, 1982.

to within numerical tolerance (currently 1e-7).

An equivalent MAD-X input file is provided in "dipedge.madx".  The user can check the matrix elements produced in "sectormap.dat" against those

obtained from ImpactX.


Run
---

This example can be run as:

* **Python** script: ``python3 run_dipedge_linear_map.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_dipedge_linear_map.py
          :language: python3
          :caption: You can copy this file from ``examples/edge_effects/run_dipedge_linear_map.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_dipedge_linear_map.py``
      
   .. literalinclude:: analysis_dipedge_linear_map.py
      :language: python3
      :caption: You can copy this file from ``examples/edge_effects/analysis_dipedge_linear_map.py``.

