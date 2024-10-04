.. _examples-fodo-channel:

FODO Channel
============

A 300m channel of 100 stable FODO cells (3m each) with a zero-current phase advance of 67.8 degrees.

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 2.82161941` m
* :math:`\alpha_\mathrm{x} = -1.59050035`
* :math:`\beta_\mathrm{y} = 2.82161941` m
* :math:`\alpha_\mathrm{y} = 1.59050035`

We use a 2 GeV electron beam with initial unnormalized rms emittance of 2 nm.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.
This test also demonstrates the ``period_sample_intervals`` capability of our beam monitor diagnostics, only creating output every 10th FODO cell


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo/run_fodo.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo/input_fodo.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo.py``

   .. literalinclude:: analysis_fodo.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/analysis_fodo.py``.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``plot_fodo.py``

   .. literalinclude:: plot_fodo.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/plot_fodo.py``.

.. figure:: https://gist.githubusercontent.com/ax3l/8ae7dcb9e07c361e002fa56d6b16cb16/raw/cc952670bb946cd7a62282bc7aa3f03f3d5faa16/fodo_channel.png
   :alt: preserved emittance in the FODO channel.

   FODO transverse emittance evolution (preserved)
