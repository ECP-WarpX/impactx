.. _examples-positron:

Positron Channel
================

Acceleration of a positron beam with large (10%) energy spread, from 10 GeV
to 2.5 TeV, for parameters based on possible staging of a laser-wakefield accelerator.

The lattice consists of 250 periods, each consisting of a quadrupole triplet
followed by 10 GeV energy gain in a uniform field.

We use a 190 pC positron beam with initial normalized rms emittance of 10 nm,
rms beam size of 5 microns, and a triangular current pulse with end-to-end pulse
length of 0.12 ps.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_positron.py`` or
* ImpactX **executable** using an input file: ``impactx input_positron.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_positron.py
          :language: python3
          :caption: You can copy this file from ``examples/positron_channel/run_positron.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_positron.in
          :language: ini
          :caption: You can copy this file from ``examples/positron_channel/input_positron.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_positron.py``

   .. literalinclude:: analysis_positron.py
      :language: python3
      :caption: You can copy this file from ``examples/positron_channel/analysis_positron.py``.
