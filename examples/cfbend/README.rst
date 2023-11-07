.. _examples-cfbend:

Combined Function Bend
======================

A single combined function bending magnet (an ideal sector bend with an upright quadrupole field component added).  The magnet
parameters are based a single CSBEND element appearing in the ELEGANT input file for the ALS-U lattice.

The beam parameters are based on:
C. Steier et al, "Status of the Conceptual Design of ALS-U", IPAC2017, WEPAB104, `DOI:10.18429/JACoW-IPAC2017-WEPAB104 <https://doi.org/10.18429/JACoW-IPAC2017-WEPAB104>`__ (2017).

A 2 GeV electron bunch with normalized transverse rms emittance of 50 pm undergoes a 3.76 deg bend.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_cfbend.py`` or
* ImpactX **executable** using an input file: ``impactx input_cfbend.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_cfbend.py
          :language: python3
          :caption: You can copy this file from ``examples/cfbend/run_cfbend.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_cfbend.in
          :language: ini
          :caption: You can copy this file from ``examples/cfbend/input_cfbend.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_cfbend.py``

   .. literalinclude:: analysis_cfbend.py
      :language: python3
      :caption: You can copy this file from ``examples/cfbend/analysis_cfbend.py``.
