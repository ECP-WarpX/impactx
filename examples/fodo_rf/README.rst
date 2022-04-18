.. _examples-fodo:

FODO Cell with RF
===================

Stable FODO cell with short RF (buncher) cavities added for
longitudinal focusing.  The phase advance in all three phase planes is
between 86-89 degrees.

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 9.80910407` m
* :math:`\alpha_\mathrm{x} = 0.0`
* :math:`\beta_\mathrm{y} = 1.31893788` m
* :math:`\alpha_\mathrm{y} = 0.0`
* :math:`\beta_\mathrm{t} = 4.6652668782` m
* :math:`\alpha_\mathrm{t} = 0.0`

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1
mm-mrad in all three phase planes.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.



.. literalinclude:: input_fodo_rf.in
   :language: ini
   :caption: You can copy this file from ``examples/fodo_rf/input_fodo_rf.in``.
