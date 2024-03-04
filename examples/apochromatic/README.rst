.. _examples-apochromat:

Apochromatic Drift-Quadrupole Beamline
======================================

Electron beam matched to the 1st-order apochromatic drift-quadrupole beamline appearing in Fig. 4a of:
C. A. Lindstrom and E. Adli, "Design of general apochromatic drift-quadrupole beam lines," Phys. Rev. Accel. Beams 19, 071002 (2016).

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 0.325` m
* :math:`\alpha_\mathrm{x} = 0`
* :math:`\beta_\mathrm{y} = 0.325` m
* :math:`\alpha_\mathrm{y} = 0`

We use a 100 GeV electron beam with an initially 6D Gaussian distribution of normalized rms emittance 1 micron and relative energy spread of 1%.

The second moments of the particle distribution after the focusing beamline should coincide with the second moments of the particle distribution before the beamline, to within the level expected due to noise due to statistical sampling.
The emittance growth due to chromatic effects remain below 1%.  In the absence of chromatic correction, the projected emittance growth is near 10%.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_apochromatic.py`` or
* ImpactX **executable** using an input file: ``impactx input_apochromatic.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_apochromatic.py
          :language: python3
          :caption: You can copy this file from ``examples/apochromatic/run_apochromatic.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_apochromatic.in
          :language: ini
          :caption: You can copy this file from ``examples/apochromatic/input_apochromatic.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_apochromatic.py``

   .. literalinclude:: analysis_apochromatic.py
      :language: python3
      :caption: You can copy this file from ``examples/apochromatic/analysis_apochromatic.py``.


.. _examples-apochromat_pl:

Apochromatic Drift-Plasma Lens Beamline
=======================================

Electron beam matched to the 3rd-order apochromatic drift-plasma lens beamline appearing in Fig. 4b of:
C. A. Lindstrom and E. Adli, "Design of general apochromatic drift-quadrupole beam lines," Phys. Rev. Accel. Beams 19, 071002 (2016).

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 0.325` m
* :math:`\alpha_\mathrm{x} = 0`
* :math:`\beta_\mathrm{y} = 0.325` m
* :math:`\alpha_\mathrm{y} = 0`

We use a 100 GeV electron beam with an initially 6D Gaussian distribution of normalized rms emittance 1 micron and relative energy spread of 1%.

The second moments of the particle distribution after the focusing beamline should coincide with the second moments of the particle distribution before the beamline, to within the level expected due to noise due to statistical sampling.
The emittance growth due to chromatic effects remain below 1%.  In the absence of chromatic correction, the projected emittance growth is near 10%.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_apochromatic_pl.py`` or
* ImpactX **executable** using an input file: ``impactx input_apochromatic_pl.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_apochromatic_pl.py
          :language: python3
          :caption: You can copy this file from ``examples/apochromatic/run_apochromatic_pl.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_apochromatic_pl.in
          :language: ini
          :caption: You can copy this file from ``examples/apochromatic/input_apochromatic_pl.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_apochromatic_pl.py``

   .. literalinclude:: analysis_apochromatic_pl.py
      :language: python3
      :caption: You can copy this file from ``examples/apochromatic/analysis_apochromatic_pl.py``.
