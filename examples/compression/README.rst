.. _examples-compression:

Ballistic Compression Using a Short RF Element
==============================================

A 20 MeV electron beam propagates through a short RF element near zero-crossing, inducing a head-tail energy correlation.
This is followed by ballistic motion in a drift, which is used to compress the rms bunch length from 16 ps to 10 ps (compression of 5/3).

The beam is not exactly on-crest (phase = -89.5 deg), so there is an energy gain of 4.5 MeV.

The transverse emittance is sufficiently small that the horizontal and verticle beam size are essentially unchanged.  Due to RF curvature, there is some growth of the longitudinal emittance.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_compression.py`` or
* ImpactX **executable** using an input file: ``impactx input_compression.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_compression.py
          :language: python3
          :caption: You can copy this file from ``examples/compression/run_compression.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_compression.in
          :language: ini
          :caption: You can copy this file from ``examples/compression/input_compression.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_compression.py``

   .. literalinclude:: analysis_compression.py
      :language: python3
      :caption: You can copy this file from ``examples/compression/analysis_compression.py``.
