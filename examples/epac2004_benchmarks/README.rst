.. _examples-fodo-rf-sc:

Cold Beam in a FODO Channel with RF Cavities (and Space Charge)
===============================================================

This example is based on the subsection of the same name in:
R. D. Ryne et al, "A Test Suite of Space-Charge Problems for Code Benchmarking", in Proc. EPAC2004, Lucerne, Switzerland.

See additional documentation in:
C. E. Mitchell et al, "ImpactX Modeling of Benchmark Tests for Space Charge Validation", in Proc. HB2023, Geneva, Switzerland.

A cold (zero momentum spread), uniform density, 250 MeV, 143 pC proton bunch propagates in a FODO lattice with 700 MHz RF
cavities added for longitudinal confinement.  The on-axis profile of the RF electric field is given by:

.. math::

   E(z)=\exp(-(4z)^4)\cos(\frac{5\pi}{2}\tanh(5z)).

The beam is matched to the 3D focusing, with space charge, using the rms envelope equations.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_rf_SC.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_rf_SC.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_rf_SC.py
          :language: python3
          :caption: You can copy this file from ``examples/epac2004_benchmarks/run_fodo_rf_SC.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_rf_SC.in
          :language: ini
          :caption: You can copy this file from ``examples/epac2004_benchmarks/input_fodo_rf_SC.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_rf_SC.py``

   .. literalinclude:: analysis_fodo_rf_SC.py
      :language: python3
      :caption: You can copy this file from ``examples/epac2004_benchmarks/analysis_fodo_rf_SC.py``.



.. _examples-thermal-beam:

Thermal Beam in a Constant Focusing Channel (with Space Charge)
===============================================================

This example is based on the subsection of the same name in:
R. D. Ryne et al, "A Test Suite of Space-Charge Problems for Code Benchmarking", in Proc. EPAC2004, Lucerne, Switzerland.

See additional documentation in:
C. E. Mitchell et al, "ImpactX Modeling of Benchmark Tests for Space Charge Validation", in Proc. HB2023, Geneva, Switzerland.

This example illustrates a stationary solution of the Vlasov-Poisson equations with spherical symmetry (in the beam
rest frame).  The distribution represents a thermal equilibrium of the form:

.. math::

   f=C\exp(-H/kT),

where :math:`C` and :math:`kT` are constants, and :math:`H` denotes the self-consistent Hamiltonian with space charge.

In this example, a 0.1 MeV, 143 pC proton bunch with :math:`kT=36\times 10^{-6}` propagates in a constant focusing lattice
with 3D isotropic focusing.  (The isotropy is exact in the beam rest frame.)

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_thermal.py`` or
* ImpactX **executable** using an input file: ``impactx input_thermal.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_thermal.py
          :language: python3
          :caption: You can copy this file from ``examples/epac2004_benchmarks/run_thermal.py``.


   .. tab-item:: Executable: Input File

       .. literalinclude:: input_thermal.in
          :language: ini
          :caption: You can copy this file from ``examples/epac2004_benchmarks/input_thermal.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_thermal.py``

   .. literalinclude:: analysis_thermal.py
      :language: python3
      :caption: You can copy this file from ``examples/epac2004_benchmarks/analysis_thermal.py``.



.. _examples-bithermal-beam:

Bithermal Beam in a Constant Focusing Channel (with Space Charge)
=================================================================

This example is based on the subsection of the same name in:
R. D. Ryne et al, "A Test Suite of Space-Charge Problems for Code Benchmarking", in Proc. EPAC2004, Lucerne, Switzerland.

See additional documentation in:
C. E. Mitchell et al, "ImpactX Modeling of Benchmark Tests for Space Charge Validation", in Proc. HB2023, Geneva, Switzerland.

This example illustrates a stationary solution of the Vlasov-Poisson equations with spherical symmetry (in the beam rest frame).
It provides a self-consistent model of a 3D bunch with a nontrivial core-halo distribution.

The distribution represents a bithermal stationary distribution of the form:

.. math::

   f=c_1\exp(-H/kT_1)+c_2\exp(-H/kT_2),

where :math:`c_j`, :math:`kT_j` :math:`(j=1,2)` are constants, and :math:`H` denotes the self-consistent Hamiltonian with space charge.

In this example, a 0.1 MeV, 143 pC proton bunch with :math:`kT_1=36\times 10^{-6}` and :math:`kT_1=900\times 10^{-6}` propagates in a constant focusing lattice
with 3D isotropic focusing.
(The isotropy is exact in the beam rest frame.)
5% of the total charge lies in the second (halo) population.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_bithermal.py`` or
* ImpactX **executable** using an input file: ``impactx input_bithermal.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_bithermal.py
          :language: python3
          :caption: You can copy this file from ``examples/epac2004_benchmarks/run_bithermal.py``.


   .. tab-item:: Executable: Input File

       .. literalinclude:: input_bithermal.in
          :language: ini
          :caption: You can copy this file from ``examples/epac2004_benchmarks/input_bithermal.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_bithermal.py``

   .. literalinclude:: analysis_bithermal.py
      :language: python3
      :caption: You can copy this file from ``examples/epac2004_benchmarks/analysis_bithermal.py``.


Visualize
---------

You can run the following script to visualize the initial and final beam distribution:

.. dropdown:: Script ``plot_bithermal.py``

   .. literalinclude:: plot_bithermal.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/plot_bithermal.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/294003440-b16185c7-2573-48d9-8998-17e116721ab5.png
   :alt: Initial and final beam distribution when running with full resolution (see inline comments in the input file/script). The bithermal distribution should stay static in this test.

   Initial and final beam distribution when running with full resolution (see inline comments in the input file/script).
   The bithermal distribution should stay static in this test.
