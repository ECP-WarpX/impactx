.. _examples-dogleg:

Dogleg
======

This is a 2-bend dogleg lattice obtained by taking the first 1/2 of the Berlin-Zeuthen magnetic bunch compression chicane:
https://www.desy.de/csr/

The primary purpose is to benchmark the reduced beam diagnostics in lattice regions with nonzero dispersion.

All parameters can be found online.  A 5 GeV electron bunch with normalized transverse rms emittance of 1 um is used.

The final expected dispersion is 267 mm.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must
agree with nominal values.

In addition, the initial and final values of :math:`\alpha_x`, :math:`\alpha_y`, :math:`\beta_x`, :math:`\beta_y`, :math:`\dispersion_x`, and :math:`\dispersion_px` must
agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_dogleg.py`` or
* ImpactX **executable** using an input file: ``impactx input_dogleg.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_dogleg.py
          :language: python3
          :caption: You can copy this file from ``examples/chicane/run_dogleg.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_dogleg.in
          :language: ini
          :caption: You can copy this file from ``examples/chicane/input_dogleg.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_dogleg.py``

   .. literalinclude:: analysis_dogleg.py
      :language: python3
      :caption: You can copy this file from ``examples/y/analysis_dogleg.py``.
