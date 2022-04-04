.. _examples-kurth:

Kurth Distribution in a Constant Focusing Channel
=================================================

Stationary Kurth distribution in a constant focusing channel (without space charge).

The distribution is radially symmetric in (x,y,t) space, and matched to a
radially symmetric constant linear focusing.

We use a 2 GeV proton beam with initial unnormalized rms emittance of 1 um
in all three phase planes.

The particle distribution should remain unchanged, to within the level expected due to numerical particle noise.
This fact is independent of the length of the channel.  This is tested using the second moments of the distribution.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


.. literalinclude:: input_kurth.in
   :language: ini
   :caption: You can copy this file from ``examples/kurth/input_kurth.in``.
