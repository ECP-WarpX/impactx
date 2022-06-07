.. _examples-kurth:

Kurth Distribution in a Constant Focusing Channel
=================================================

Stationary Kurth distribution in a constant focusing channel (without space charge).

The distribution is radially symmetric in (x,y,t) space, and matched to a
radially symmetric constant linear focusing.

We use a 2 GeV proton beam with initial unnormalized rms emittance of 1 um
in all three phase planes.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This fact is independent of the length of the channel.  This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as a Python script (``python3 run_kurth.py``) or with an app with an input file (``impactx input_kurth.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_kurth.py
          :language: python3
          :caption: You can copy this file from ``examples/kurth/run_kurth.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_kurth.in
          :language: ini
          :caption: You can copy this file from ``examples/kurth/input_kurth.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_kurth.py``

   .. literalinclude:: analysis_kurth.py
      :language: python3
      :caption: You can copy this file from ``examples/kurth/analysis_kurth.py``.


.. _examples-kurth-10nC:

Kurth Distribution in a Constant Focusing Channel with Space Charge
===================================================================

Stationary Kurth distribution in a constant focusing channel with space charge.

The distribution is radially symmetric in (x,y,t) space, and matched to a
radially symmetric constant linear focusing.

We use a 2 GeV proton beam with initial unnormalized rms emittance of 1 um
in all three phase planes.  The bunch charge is set to 10 nC, corresponding
to a transverse tune depression ratio of 0.67.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This fact is independent of the length of the channel.  This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`


Run
---

This example can be run as an app with an input file (``impactx input_kurth_periodic.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script (TODO)

       This one is to-do.

   .. tab-item:: App Input File

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
