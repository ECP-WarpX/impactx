.. _examples-fodo-tune:

Tune Calculation in a Periodic FODO Channel
===========================================

This is identical to the :ref:`FODO example <examples-fodo>`, except that tracking for 100 periods is used to extract the horizontal tune.

Stable FODO cell with a zero-current phase advance of 67.8 degrees, corresponding to Qx = Qy = 0.1883.

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 2.82161941` m
* :math:`\alpha_\mathrm{x} = -1.59050035`
* :math:`\beta_\mathrm{y} = 2.82161941` m
* :math:`\alpha_\mathrm{y} = 1.59050035`

We use a 2 GeV electron beam with initial unnormalized rms emittance of 2 nm.

The horizontal tune of a single particle is obtained from period-by-period tracking data using several algorithms.

In this test, the computed horizontal tune must agree with the nominal value to within acceptable tolerance.

This example requires installation of `PyNAFF <https://github.com/nkarast/PyNAFF>`__.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_tune.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_tune.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_tune.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo_tune/run_fodo_tune.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_tune.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo_tune/input_fodo_tune.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_tune.py``

   .. literalinclude:: analysis_fodo_tune.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_tune/analysis_fodo_tune.py``.
