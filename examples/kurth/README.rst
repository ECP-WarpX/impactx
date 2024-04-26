.. _examples-kurth-periodic:

Kurth Distribution in a Periodic Focusing Channel
=================================================

Matched Kurth distribution in a periodic focusing channel (without space charge).

The distribution is radially symmetric in (x,y,t) space, and matched to a
radially symmetric periodic linear focusing lattice with a phase advance of 121 degrees.

We use a 2 GeV proton beam with initial unnormalized rms emittance of 1 um
in all three phase planes.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_kurth_periodic.py`` or
* ImpactX **executable** using an input file: ``impactx input_kurth_periodic.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_kurth_periodic.py
          :language: python3
          :caption: You can copy this file from ``examples/kurth/run_kurth_periodic.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_kurth_periodic.in
          :language: ini
          :caption: You can copy this file from ``examples/kurth/input_kurth_periodic.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_kurth_periodic.py``

   .. literalinclude:: analysis_kurth_periodic.py
      :language: python3
      :caption: You can copy this file from ``examples/kurth/analysis_kurth_periodic.py``.


.. _examples-kurth-10nC-periodic:

Kurth Distribution in a Periodic Focusing Channel with Space Charge
===================================================================

Matched Kurth distribution in a periodic focusing channel with space charge.

The distribution is radially symmetric in (x,y,t) space, and matched to a
radially symmetric constant linear focusing.

We use a 2 GeV proton beam with initial unnormalized rms emittance of 1 um
in all three phase planes.  The bunch charge is set to 10 nC, depressing the
phase advance from 121 degrees to 74 degrees.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`


Run
---

This example can be run as a Python script (``python3 run_kurth_10nC_periodic.py``) or as an app with an input file (``impactx input_kurth_10nC_periodic.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_kurth_10nC_periodic.py
          :language: python3
          :caption: You can copy this file from ``examples/kurth/run_kurth_10nC_periodic.py``.


   .. tab-item:: App Input File

       .. literalinclude:: input_kurth_10nC_periodic.in
          :language: ini
          :caption: You can copy this file from ``examples/kurth/input_kurth_10nC_periodic.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_kurth_10nC_periodic.py``

   .. literalinclude:: analysis_kurth_10nC_periodic.py
      :language: python3
      :caption: You can copy this file from ``examples/kurth/analysis_kurth_10nC_periodic.py``.
