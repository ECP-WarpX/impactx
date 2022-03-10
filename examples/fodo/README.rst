.. _examples-fodo:

FODO Cell
=========

Stable FODO cell with a zero-current phase advance of 67.8 degrees.

The matched Twiss parameters at entry are:

horizontal beta = 2.82161941 m
horizontal alpha = -1.59050035

vertical beta = 2.82161941 m
vertical alpha = 1.59050035

We use a 2 GeV proton beam with initial unnormalized rms emittance of 2 nm.

The second moments of the particle distribution after the FODO cell
should coincide with the second moments of the particle distribution
before the FODO cell, to within the level expected due to noise due
to statistical sampling.

In this test, the initial and final values of sigma_x, sigma_y, sigma_t,
emittance_x, emittance_y, and emittance_t must agree with nominal values.


.. literalinclude:: input_fodo.in
   :language: ini
   :caption: You can copy this file from ``examples/fodo/input_fodo.in``.
