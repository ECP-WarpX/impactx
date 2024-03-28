.. _examples-triplet:

Optimized Triplet
=================

Optimization of focusing parameters for a quadrupole triplet.
A 2 GeV electron beam is strongly focused from lattice initial parameters:

.. math::

   \beta_x = \beta_y = 40\,\mathrm{m}\\
   \alpha_x = -\alpha_y = 2

to final lattice parameters:

.. math::

   \beta_x = \beta_y = 0.55\,\mathrm{m}\\
   \alpha_x = \alpha_y = 0

resulting in a reduction by a factor of 8.5 in the horizontal and vertical beam sizes.

Here, we start with a desired spatial layout of the triplet and find the quadrupole strengths through numerical optimization (by minimizing the L2 norm of alpha and beta) over multiple ImpactX simulations.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

.. tab-set::

   .. tab-item:: SciPy Optimizers

      This example uses `scipy.optimize <https://docs.scipy.org/doc/scipy/reference/optimize.html>`__ (methods: `Nelder-Mead <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html>`__ or `L-BFGS-B <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html>`__) to find the quadrupole strengths by minimizing the objective.
      Conventional optimization algorithms like this work best if there is only a global minima in the objective.

      This example can be run as a **Python** script: ``python3 run_triplet.py``.

      .. literalinclude:: run_triplet.py
         :language: python3
         :caption: You can copy this file from ``examples/optimize_triplet/run_triplet.py``.

   .. tab-item:: Xopt

      This example uses `Xopt <https://christophermayes.github.io/Xopt>`__ (methods: `Nelder-Mead <https://christophermayes.github.io/Xopt/examples/scipy/neldermead>`__ or `TuRBO <https://christophermayes.github.io/Xopt/examples/single_objective_bayes_opt/turbo_tutorial>`__) to find the quadrupole strengths by minimizing the objective.

      Conventional optimization algorithms like ``Nelder-Mead`` work best if there is only a global minima in the objective.
      Machine-learning based, surrogate optimization works well for highly dimensional inputs and/or to find global minima in an objective that has potentially many local minima, where conventional optimizers can get stuck.
      At the same time, the ML method Bayesian Optimization (BO) is prone to over-explore an objective (at the cost of finding a point closer to the global minima).
      The variation ``Trust Region Bayesian Optimization (TuRBO)`` was developed to narrow down further on found minima.

      This example can be run as a **Python** script: ``python3 tests/python/test_xopt.py``.

      .. literalinclude:: test_xopt.py
         :language: python3
         :caption: You can copy this file from ``tests/python/test_xopt.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_triplet.py``

   .. literalinclude:: analysis_triplet.py
      :language: python3
      :caption: You can copy this file from ``examples/optimize_triplet/analysis_triplet.py``.


Visualize
---------

You can run the following script to visualize the optimized beam evolution over time:

.. dropdown:: Script ``plot_triplet.py``

   .. literalinclude:: plot_triplet.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/plot_triplet.py``.

.. figure:: https://gist.github.com/assets/1353258/d3bf0184-a102-47cf-98d3-da2a5d511b96
   :alt: CS Twiss beta

   CS Twiss beta
