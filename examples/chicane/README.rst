.. _examples-chicane:

Chicane
=======

Berlin-Zeuthen magnetic bunch compression chicane:
https://www.desy.de/csr/

All parameters can be found online.
A 5 GeV electron bunch with normalized transverse rms emittance of 1 um undergoes longitudinal compression by a factor of 10 in a standard 4-bend chicane.

The emittances should be preserved, and the rms pulse length should decrease by the compression factor (10).

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example an either be run with an inputs file (``impactx input_chicane.in``):

.. literalinclude:: input_chicane.in
   :language: ini
   :caption: You can copy this file from ``examples/chicane/input_chicane.in``.


Analyze
-------

We run the following script to analyze correctness:

.. literalinclude:: analysis_chicane.py
   :language: python3
   :caption: You can copy this file from ``examples/chicane/analysis_chicane.py``.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. literalinclude:: plot_chicane.py
   :language: python3
   :caption: You can copy this file from ``examples/chicane/plot_chicane.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/180332191-f9ce11fc-8c56-4713-a91a-2ad12ab09805.png
 :alt: Chicane beam width and emittance evolution

.. figure:: https://user-images.githubusercontent.com/1353258/180332189-a0e6933b-d3c6-4170-acd0-d79b366602e0.png
 :alt: Chicane phase space evolution
