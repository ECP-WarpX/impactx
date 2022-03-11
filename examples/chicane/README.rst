.. _examples-chicane:

Chicane
=======

Berlin-Zeuthen magnetic bunch compression chicane:
https://www.desy.de/csr/

All parameters can be found online.
A 5 GeV electron bunch with normalized transverse rms emittance of 1 um undergoes longitudinal compression by a factor of 10 in a standard 4-bend chicane.

The emittances should be preserved, and the rms pulse length should decrease by the compression factor (10).

In this test, the initial and final values of sigma_x, sigma_y, sigma_t,
emittance_x, emittance_y, and emittance_t must agree with nominal values.


.. literalinclude:: input_chicane.in
   :language: ini
   :caption: You can copy this file from ``examples/chicane/input_chicane.in``.
