.. _examples-expanding:

Expanding Beam in Free Space
============================

A coasting bunch expanding freely in free space under its own space charge.

We use a cold (zero emittance) 250 MeV electron bunch whose
initial distribution is a uniformly-populated 3D ball of radius R0 = 1 mm when viewed in the bunch rest
frame.

In the laboratory frame, the bunch is a uniformly-populated ellipsoid, which
expands to twice its original size.  This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

This test uses mesh-refinement to solve the space charge force.
The coarse grid wraps the beam maximum extent by 300%, emulating "open boundary" conditions.
The refined grid in level 1 spans 110% of the beam maximum extent.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_expanding.py`` or
* ImpactX **executable** using an input file: ``impactx input_expanding.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. literalinclude:: run_expanding.py
         :language: python3
         :caption: You can copy this file from ``examples/expanding/run_expanding.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_expanding.in
          :language: ini
          :caption: You can copy this file from ``examples/expanding/input_expanding.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_expanding.py``

   .. literalinclude:: analysis_expanding.py
      :language: python3
      :caption: You can copy this file from ``examples/expanding/analysis_expanding.py``.
