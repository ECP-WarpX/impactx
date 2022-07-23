.. _examples-fodo:

FODO Cell
=========

Stable FODO cell with a zero-current phase advance of 67.8 degrees.

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 2.82161941` m
* :math:`\alpha_\mathrm{x} = -1.59050035`
* :math:`\beta_\mathrm{y} = 2.82161941` m
* :math:`\alpha_\mathrm{y} = 1.59050035`

We use a 2 GeV proton beam with initial unnormalized rms emittance of 2 nm.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example an either be run with an inputs file (``impactx input_fodo.in``):

.. literalinclude:: input_fodo.in
   :language: ini
   :caption: You can copy this file from ``examples/fodo/input_fodo.in``.

Or as a Python script (``python3 run_fodo.py``):

.. literalinclude:: run_fodo.py
   :language: python3
   :caption: You can copy this file from ``examples/fodo/run_fodo.py``.

Both execution modes can also be prefixed with an MPI executor, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.


Analyze
-------

We run the following script to analyze correctness:

.. literalinclude:: analysis_fodo.py
   :language: python3
   :caption: You can copy this file from ``examples/fodo/analysis_fodo.py``.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. literalinclude:: plot_fodo.py
   :language: python3
   :caption: You can copy this file from ``examples/fodo/plot_fodo.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/180287840-8561f6fd-278f-4856-abd8-04fbdb78c8ff.png
 :alt: FODO beam width and emittance evolution

.. figure:: https://user-images.githubusercontent.com/1353258/180287845-eb0210a7-2500-4aa9-844c-67fb094329d3.png
 :alt: FODO phase space evolution
