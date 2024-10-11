.. _examples-aperture:

Aperture Collimation
====================

Proton beam undergoing collimation by a rectangular boundary aperture.


We use a 250 MeV proton beam with a horizontal rms beam size of 1.56 mm and a vertical rms beam size of 2.21 mm.

After a short drift of 0.123, the beam is scraped by a 1 mm x 1.5 mm rectangular aperture.

In this test, the initial values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.
The test fails if:

* any of the final coordinates for the valid (not lost) particles lie outside the aperture boundary or
* any of the lost particles are inside the aperture boundary or
* if the sum of lost and kept particles is not equal to the initial particles or
* if the recorded position :math:`s` for the lost particles does not coincide with the drift distance.


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


.. _examples-aperture-pepperpot:

Aperture Collimation with Periodic Masking
===========================================

Proton beam undergoing collimation by a periodic array of rectangular apertures, such as those used in a pepperpot emittance measurement.

We use a 250 MeV proton beam with a horizontal rms beam size of 1.56 mm and a vertical rms beam size of 2.21 mm.

After a short drift of 0.123 m, the beam intercepts a periodic array of apertures of period 1 mm (in the horizontal and vertical), each of which is 0.15 mm x 0.1 mm in size.

In this test, the initial values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nomin>
The test fails if:

* any of the final coordinates for the valid (not lost) particles lie outside the aperture boundary or
* any of the lost particles are inside the aperture boundary or
* if the sum of lost and kept particles is not equal to the initial particles


Run
---

This example can be run as a Python script (``python3 run_aperture_pepperpot.py``) or with an app with an input file (``impactx input_aperture_pepperpot.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_aperture_pepperpot.py
          :language: python3
          :caption: You can copy this file from ``examples/aperture/run_aperture_pepperpot.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_aperture_pepperpot.in
          :language: ini
          :caption: You can copy this file from ``examples/aperture/input_aperture_pepperpot.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_aperture_pepperpot.py``

   .. literalinclude:: analysis_aperture_pepperpot.py
      :language: python3
      :caption: You can copy this file from ``examples/aperture/analysis_aperture_pepperpot.py``.
