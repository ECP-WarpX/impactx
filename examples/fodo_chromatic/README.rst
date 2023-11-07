.. _examples-fodo-chr:

FODO Cell, Chromatic
====================

Stable FODO cell with a zero-current phase advance of 67.8 degrees,
with chromatic focusing effects included.

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 2.82161941` m
* :math:`\alpha_\mathrm{x} = -1.59050035`
* :math:`\beta_\mathrm{y} = 2.82161941` m
* :math:`\alpha_\mathrm{y} = 1.59050035`

We use a 2 GeV electron beam with initial unnormalized rms emittance of 2 nm.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_chr.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_chr.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_chr.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo/run_fodo_chr.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_chr.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo/input_fodo_chr.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_chr.py``

   .. literalinclude:: analysis_fodo_chr.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/analysis_fodo_chr.py``.
