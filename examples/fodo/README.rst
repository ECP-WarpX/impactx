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



.. literalinclude:: input_fodo.in
   :language: ini
   :caption: You can copy this file from ``examples/fodo/input_fodo.in``.
