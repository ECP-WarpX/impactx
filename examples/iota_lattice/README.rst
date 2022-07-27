.. _examples-iotalattice:

The "bare" linear lattice of the Fermilab IOTA storage ring.
=============================================================

The linear lattice of the IOTA storage ring, configured for operation with a
2.5 MeV proton beam.

The drift regions available for insertion of the special nonlinear magnetic 
element for integrable optics experiments are denoted "dnll".

The second moments of the particle distribution after a single turn should
coincide with the initial secton moments of the particle distribution, to
within the level expected due to numerical particle noise.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\eepsilon_t`
must agree with nominal values.


.. literalinclude:: input_iotalattice.in
   :language: ini
   :caption: You can copy this file from ``examples/iota_lattice/input_iotalattice.in``.
