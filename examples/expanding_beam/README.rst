
.. _examples-expanding:

Expanding Beam in Free Space with 3D Space Charge
=================================================

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
The grid spacing is adaptively adjusted as the beam evolves.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_expanding_fft.py`` or
* ImpactX **executable** using an input file: ``impactx input_expanding_fft.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

We also provide the same example with the multi-grid (MLMG) Poisson solver,
and with the MLMG solver inside a grounded beam pipe (:pp:param:`eb.shape`) that is far enough from the beam to not change the result.

.. tab-set::

   .. tab-item:: Python: Script (FFT)

      .. literalinclude:: run_expanding_fft.py
         :language: python3
         :caption: You can copy this file from ``examples/expanding/run_expanding_fft.py``.

   .. tab-item:: Python: Script (MLMG)

      .. literalinclude:: run_expanding_mlmg.py
         :language: python3
         :caption: You can copy this file from ``examples/expanding/run_expanding_mlmg.py``.

   .. tab-item:: Executable: Input File (FFT)

       .. literalinclude:: input_expanding_fft.in
          :language: ini
          :caption: You can copy this file from ``examples/expanding/input_expanding_fft.in``.

   .. tab-item:: Executable: Input File (MLMG)

       .. literalinclude:: input_expanding_mlmg.in
          :language: ini
          :caption: You can copy this file from ``examples/expanding/input_expanding_mlmg.in``.

   .. tab-item:: Executable: Input File (MLMG, beam pipe)

       .. literalinclude:: input_expanding_mlmg_eb.in
          :language: ini
          :caption: You can copy this file from ``examples/expanding/input_expanding_mlmg_eb.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_expanding.py``

   .. literalinclude:: analysis_expanding.py
      :language: python3
      :caption: You can copy this file from ``examples/expanding/analysis_expanding.py``.



.. _examples-expanding-fft-2d:

Expanding Beam in Free Space with 2D Space Charge
=================================================

A long, coasting unbunched beam expanding freely in free space under its own 2D space charge.

We use a cold (zero emittance) 250 MeV electron bunch whose
initial distribution is a uniformly-populated cylinder of radius R0 = 1 mm.

In the laboratory frame, the beam expands to twice its original transverse size.  This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_expanding_fft_2D.py`` or
* ImpactX **executable** using an input file: ``impactx input_expanding_fft_2D.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

We also provide the same example with the multi-grid (MLMG) Poisson solver.

.. tab-set::

   .. tab-item:: Python: Script (FFT)

      .. literalinclude:: run_expanding_fft_2D.py
         :language: python3
         :caption: You can copy this file from ``examples/expanding/run_expanding_fft_2D.py``.

   .. tab-item:: Executable: Input File (FFT)

       .. literalinclude:: input_expanding_fft_2D.in
          :language: ini
          :caption: You can copy this file from ``examples/expanding/input_expanding_fft_2D.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_expanding_fft_2D.py``

   .. literalinclude:: analysis_expanding_fft_2D.py
      :language: python3
      :caption: You can copy this file from ``examples/expanding/analysis_expanding_fft_2D.py``.



.. _examples-expanding-fft-2d-test-particles:

Expanding Beam in Free Space with 2D Space Charge, with Test Particles
======================================================================

A long, coasting unbunched beam expanding freely in free space under its own 2D space charge.

The problem parameters are identical to :ref:`the above example <examples-expanding-fft-2d>`, except that a set of test particles is loaded and tracked, together with the main particle bunch.

The test particle data is collected once per space charge slice, using the callback hook feature, and the test particle orbits are plotted vs. distance s.

Run
---

This example can be run as:

* **Python** script: ``python3 run_expanding_fft_2D_test_particles.py`` or

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script (FFT)

      .. literalinclude:: run_expanding_fft_2D_test_particles.py
         :language: python3
         :caption: You can copy this file from ``examples/expanding/run_expanding_fft_2D_test_particles.py``.

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.


.. _examples-expanding-acc:

Expanding Beam in a Uniform Accelerating Field with 3D Space Charge
===================================================================

A bunch in a uniform longitudinal electric field, expanding under its own space charge.

The problem is identical to :ref:`the above example <examples-expanding-fft>`, except for the presence of the accelerating field.
Since all forces are linear in this case, the problem is characterized by the rms envelope equations.

In this test, the final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` obtained using
particle tracking must agree with the corresponding values obtained using envelope tracking.


Run
---

This example can be run as:

* **Python** script: ``python3 run_expanding_acc_fft.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script (FFT)

      .. literalinclude:: run_expanding_acc_fft.py
         :language: python3
         :caption: You can copy this file from ``examples/expanding/run_expanding_acc_fft.py``.


Analyze
-------

Analysis is performed in-situ for this problem




.. _examples-expanding-2d-acc:

Expanding Beam in Uniform Accelerating Field with 2D Space Charge
=================================================================

A long, unbunched beam in a uniform longitudinal electric field, expanding under its own 2D space charge.

The problem is nearly identical to :ref:`the above example <examples-expanding-fft-2d>`, except for the presence of the accelerating field.
The propagation distance has also been modified.
Since all forces are linear in this case, the problem is characterized by the rms envelope equations.

In this test, the final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` obtained using
particle tracking must agree with the corresponding values obtained using envelope tracking.

Run
---

This example can be run as:

* **Python** script: ``python3 run_expanding_acc_fft_2D.py`` or

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script (FFT)

      .. literalinclude:: run_expanding_acc_fft_2D.py
         :language: python3
         :caption: You can copy this file from ``examples/expanding/run_expanding_acc_fft_2D.py``.


Analyze
-------

Analysis is performed in-situ for this problem.
