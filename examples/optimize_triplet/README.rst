.. _examples-triplet:

Optimized Triplet
=================

Optimization of focusing parameters for a quadrupole triplet.
A 2 GeV electron beam is strongly focused from lattice initial parameters:

:math:`\beta_x = \beta_y = 40` m
:math:`\alpha_x = -\alpha_y = 2`

to final lattice parameters:

:math:`\beta_x = \beta_y = 0.55` m
:math:`\alpha_x = \alpha_y = 0`

resulting in a reduction by a factor of 8.5 in the horizontal and vertical beam sizes.

Here, we start with a desired spatial layout of the triplet and find the quadrupole strengths through numerical optimization (by `minimizing the L2 norm of alpha and beta <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`__) over multiple ImpactX simulations.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as a **Python** script: ``python3 run_triplet.py``.

.. literalinclude:: run_triplet.py
  :language: python3
  :caption: You can copy this file from ``examples/optimize_triplet/run_triplet.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_triplet.py``

   .. literalinclude:: analysis_triplet.py
      :language: python3
      :caption: You can copy this file from ``examples/optimize_triplet/analysis_triplet.py``.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``plot_triplet.py``

   .. literalinclude:: plot_triplet.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/plot_triplet.py``.

.. figure:: https://private-user-images.githubusercontent.com/1353258/309984511-ba3ca6a0-e73e-4c69-9e7e-b60685433a6b.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MDk2MTM4NjQsIm5iZiI6MTcwOTYxMzU2NCwicGF0aCI6Ii8xMzUzMjU4LzMwOTk4NDUxMS1iYTNjYTZhMC1lNzNlLTRjNjktOWU3ZS1iNjA2ODU0MzNhNmIucG5nP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQVZDT0RZTFNBNTNQUUs0WkElMkYyMDI0MDMwNSUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyNDAzMDVUMDQzOTI0WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9ODhkZWMyNjkxYjQwZjgzNGNiYzU5N2UxMmM0ZGY5YTRmYTkzZGJkNmY5Mzg4NzAzMTFlMTEyMGMyNGFjZjU1OCZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QmYWN0b3JfaWQ9MCZrZXlfaWQ9MCZyZXBvX2lkPTAifQ.yfNmqJohzemrw5LQum6zeqfGPxhAi9VEBXH_ieEALro
   :alt: CS Twiss beta

   CS Twiss beta
