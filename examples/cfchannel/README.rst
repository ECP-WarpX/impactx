.. _examples-cfchannel:

Constant Focusing Channel
==========================

Stationary beam in a constant focusing channel (without space charge).

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 1.0` m
* :math:`\alpha_\mathrm{x} = 0.0`
* :math:`\beta_\mathrm{y} = 1.0` m
* :math:`\alpha_\mathrm{y} = 0.0`

We use a 2 GeV proton beam with initial unnormalized rms emittance of 1 um. 
The longitudinal beam parameters are chosen so that the bunch has radial
symmetry when viewed in the beam rest frame.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This fact is independent of the length of the channel.  This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


.. literalinclude:: input_cfchannel.in
   :language: ini
   :caption: You can copy this file from
``examples/cfchannel/input_cfchannel.in``.
