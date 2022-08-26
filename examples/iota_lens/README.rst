.. _examples-iotalens:

A nonlinear focusing channel based on the IOTA nonlinear lens
=============================================================

A constant focusing channel with nonlinear focusing, using a string of thin
IOTA nonlinear lens elements alternating with constant focusing elements.

We use a 2.5 MeV proton beam, corresponding to the nominal IOTA proton energy.

The two functions H (Hamiltonian) and I (the second invariant) should remain unchanged for all particles.

In this test, the initial and final values of :math:`\mu_H`, :math:`\sigma_H`, :math:`\mu_I`, :math:`\sigma_I` must agree with nominal values.


.. literalinclude:: input_iotalens.in
   :language: ini
   :caption: You can copy this file from ``examples/iota_lens/input_iotalens.in``.

This test can also be run as a Python script (``python3 run_iotalens.py``):

.. literalinclude:: run_iotalens.py
   :language: python3
   :caption: You can copy this file from ``examples/iota_lens/run_iotalens.py``.
