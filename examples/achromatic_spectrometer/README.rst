.. _examples-achromatic-spectrometer:

Achromatic Spectrometer
=======================

A spectrometer beamline using a bending dipole.
A :py:class:`transversely-tapered plasma lens <impactx.elements.impactx.elements.TaperedPL>` is added for chromatic correction.
The tapered plasma lens design is based on:

C. A. Lindstrom, presentation at the EuroNNAc Special Topics Workshop 2022, `slides <https://agenda.infn.it/event/28376/contributions/178724/attachments/96899/133588/Lindstr%C3%B8m,%20EuroNNAc%20workshop,%2022%20Sep%202022.pdf>`__
"Solutions and challenges for a multi-stage plasma accelerator",

https://agenda.infn.it/event/28376/contributions/178724/attachments/96899/133588/Lindstr%C3%B8m,%20EuroNNAc%20workshop,%2022%20Sep%202022.pdf

with a transverse (horizontal) taper

.. math::

   B_x = g \left( y + \frac{xy}{D_x} \right), \quad \quad B_y = -g \left(x + \frac{x^2 + y^2}{2 D_x} \right)

where :math:`g` is the (linear) field gradient in T/m and :math:`D_x` is the targeted horizontal dispersion in m.

We use a 1 GeV electron beam with initial normalized rms emittance of 2 microns and 2% rms relative energy spread.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_spectrometer.py`` or
* ImpactX **executable** using an input file: ``impactx input_spectrometer.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_spectrometer.py
          :language: python3
          :caption: You can copy this file from ``examples/achromatic_spectrometer/run_spectrometer.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_spectrometer.in
          :language: ini
          :caption: You can copy this file from ``examples/achromatic_spectrometer/input_spectrometer.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_spectrometer.py``

   .. literalinclude:: analysis_spectrometer.py
      :language: python3
      :caption: You can copy this file from ``examples/achromatic_spectrometer/analysis_spectrometer.py``.
