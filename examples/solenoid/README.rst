.. _examples-solenoid:

Solenoid channel
================

Proton beam undergoing 90 deg X-Y rotation in an ideal solenoid channel.

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 2.4321374875` m
* :math:`\alpha_\mathrm{x} = 0.0`
* :math:`\beta_\mathrm{y} = 2.4321374875` m
* :math:`\alpha_\mathrm{y} = 0.0`

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1 micron
in the horizontal plane, and 2 micron in the vertical plane.

The solenoid magnetic field corresponds to B = 2 T.

The second moments of the particle distribution after the solenoid channel are rotated by 90 degrees:  the final horizontal moments should coincide with the
initial vertical moments, and vice-versa, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_solenoid.py`` or
* ImpactX **executable** using an input file: ``impactx input_solenoid.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_solenoid.py
          :language: python3
          :caption: You can copy this file from ``examples/solenoid/run_solenoid.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_solenoid.in
          :language: ini
          :caption: You can copy this file from ``examples/solenoid/input_solenoid.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_solenoid.py``

   .. literalinclude:: analysis_solenoid.py
      :language: python3
      :caption: You can copy this file from ``examples/solenoid/analysis_solenoid.py``.
