.. _theory:

Introduction
============

.. note::

   TODO


Assumptions
-----------

This is a work-in-progress list of physical assumptions implemented in the numerics of ImpactX.


Tracking and Lattice Optics
"""""""""""""""""""""""""""

* **tracking through lattice optics:** is treated through linear order with respect to the reference particle

  * **velocity spread:** the above linearization implies that, when solving space-charge effects, we assume that the relative spread of velocities of particles in the beam is negligible compared to the velocity of the reference particle


Space Charge (Poisson Solver)
"""""""""""""""""""""""""""""

* **electrostatic in the bunch frame:** we assume there are no retardation effects and we solve the Poisson equation in the bunch frame
