
.. _examples-csr-kick:

Momentum Kick in the 1D Steady-State CSR Model
===============================================

The test compares the energy change due to CSR experienced by particles in a bending dipole within a long, initially cold, Gaussian bunch against the analytical model in:

J. Qiang, Phys. Rev. Accel. Beams 28, 114602 (2025), eqs. (31-32).

We use a cold (zero emittance) 100 MeV electron bunch in a drift space, with a single slice.

Here :math:`\sigma_x = \sigma_y`, so the predicted momentum kick can be evaluated in closed form.

In this test, the final momentum of every particle is compared against its predicted value.  Both the maximum error and its rms value over the bunch population are returned.  Both values must be zero within a small tolerance.


Run
---

This example can be run as:

* **Python** script: ``python3 run_csr_kick.py`` or ``python3 run_csr_kick.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Gauss2p5D space charge solver

      .. literalinclude:: run_sc_kick_Gauss2p5D.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/run_sc_kick_Gauss2p5D.py``.

   .. tab-item:: 2p5D (PIC) space charge solver

      .. literalinclude:: run_sc_kick_PIC2p5D.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/run_sc_kick_PIC2p5D.py``.

Analyze
-------

We run the following script to analyze correctness:

.. tab-set::

   .. tab-item:: Gauss2p5D space charge solver

      .. literalinclude:: analysis_sc_kick_Gauss2p5D.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/analysis_sc_kick_Gauss2p5D.py``.

   .. tab-item:: 2p5D (PIC) space charge solver

      .. literalinclude:: run_sc_kick_PIC2p5D.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/analysis_sc_kick_PIC2p5D.py``.


.. _examples-kick-Gauss2p5D-long-kick-off:

Turning Off the Longitudinal Kick in the 2.5D Gaussian Space Charge Model
=========================================================================

This test is identical to the test examples-kick-Gauss2p5D, except that the longitudinal space charge kick is turned off.

The initial and final values of pt should agree to within a small tolerance.



.. _examples-kick-Gauss2p5D-protons:

Momentum Kick in the 2.5D Gaussian Space Charge Model (Protons)
===============================================================

This test is based on the test :ref:`examples-kick-Gauss2p5D`, except that the electrons bunch is replaced by a proton bunch, resulting in a beam with low relativistic beta.

The purpose is to validate the ``Gauss2p5D`` model for low-energy proton beams.

The test compares the space charge momentum kick experienced by particles within a long, initially cold, Gaussian bunch against the 2.5D analytical model in:

`J. Qiang, Phys. Rev. Accel. Beams 28, 114602 (2025), eqs. (31-32). <https://doi.org/10.1103/npx1-nclt>`__

We use a cold (zero emittance) 100 MeV proton bunch in a drift space, with a single space charge slice.

Here :math:`\sigma_x = \sigma_y`, so the predicted momentum kick can be evaluated in closed form.

In this test, the final momentum of every particle is compared against its predicted value.  Both the maximum error and its rms value over the bunch population are returned.  Both values must be zero within a small tolerance.


Run
---

This example can be run as:

* **Python** script: ``python3 run_sc_kick_Gauss2p5D_protons.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Gauss2p5D space charge solver

      .. literalinclude:: run_sc_kick_Gauss2p5D_protons.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/run_sc_kick_Gauss2p5D_protons.py``.

Analyze
-------

We run the following script to analyze correctness:

.. tab-set::

   .. tab-item:: Gauss2p5D space charge solver

      .. literalinclude:: analysis_sc_kick_Gauss2p5D_protons.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/analysis_sc_kick_Gauss2p5D_protons.py``.

Visualize
---------

You can run the following script to visualize the longitudinal kick applied to the particles in the bunch:

.. dropdown:: Script ``plot_sc_kick_Gauss2p5D_protons.py``

   .. literalinclude:: plot_sc_kick_Gauss2p5D_protons.py
      :language: python3
      :caption: You can copy this file from ``examples/space_charge_gaussian_kick/plot_sc_kick_Gauss2p5D_protons.py``.



.. _examples-kick-PIC2p5D-protons:

Momentum Kick in the 2.5D PIC Space Charge Model (Protons)
==========================================================

This test is based on the test :ref:`examples-kick-Gauss2p5D-protons`, except that the 2.5D FFT-based particle-in-cell solver is used.

The purpose is to validate the ``2p5D`` model for low-energy proton beams.

The test compares the space charge momentum kick experienced by particles within a long, initially cold, Gaussian bunch against the 2.5D analytical model in:

`J. Qiang, Phys. Rev. Accel. Beams 28, 114602 (2025), eqs. (31-32). <https://doi.org/10.1103/npx1-nclt>`__

We use a cold (zero emittance) 100 MeV electron bunch in a drift space, with a single space charge slice.

Here :math:`\sigma_x = \sigma_y`, so the predicted momentum kick can be evaluated in closed form.

In this test, the final momentum of every particle is compared against its predicted value.  Both the maximum error and its rms value over the bunch population are returned.  Both values must be zero within a small tolerance.


Run
---

This example can be run as:

* **Python** script: ``python3 run_sc_kick_PIC2p5D_protons.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: 2p5D (PIC) space charge solver

      .. literalinclude:: run_sc_kick_PIC2p5D_protons.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/run_sc_kick_PIC2p5D_protons.py``.

Analyze
-------

We run the following script to analyze correctness:

.. tab-set::


   .. tab-item:: 2p5D (PIC) space charge solver

      .. literalinclude:: run_sc_kick_PIC2p5D_protons.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/analysis_sc_kick_PIC2p5D_protons.py``.
