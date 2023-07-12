.. _examples-aperture:

Aperture collimation
=====================

Proton beam undergoing collimation by a rectangular boundary aperture.


We use a 250 MeV proton beam with a horizontal rms beam size of 1.56 mm and a vertical rms beam size of 2.21 mm.

The beam is scraped by a 1 mm x 1.5 mm rectangular aperture.

In this test, the initial values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.
The test fails if any of the final coordinates for the valid (not lost) particles lie outside the aperture boundary.


Run
---

This example can be run as a Python script (``python3 run_aperture.py``) or with an app with an input file (``impactx input_aperture.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_aperture.py
          :language: python3
          :caption: You can copy this file from ``examples/aperture/run_aperture.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_aperture.in
          :language: ini
          :caption: You can copy this file from ``examples/aperture/input_aperture.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_aperture.py``

   .. literalinclude:: analysis_aperture.py
      :language: python3
      :caption: You can copy this file from ``examples/aperture/analysis_aperture.py``.
