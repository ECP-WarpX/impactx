.. _theory:

Introduction
============

Concepts
--------

Reference Trajectory
""""""""""""""""""""

ImpactX is an *s*-based beam dynamics code, evolving particles *relative* to a reference trajectory.

A reference trajectory is favorable instead of, e.g., differences to the beam centroid, because:

#. the reference trajectory is the ideal single-particle orbit used as part of the optics design, so it is better that it can be computed independently of the beam distribution,
#. differences between the beam centroid and the reference particle are important when investigating the effects of misalignments and errors for beamline designs,
#. the fields are usually specified relative to the reference trajectory, and the dynamics becomes more nonlinear as one moves away from the reference trajectory, so knowing how far the beam particles are form the reference trajectory is important for accuracy,
#. we want to keep track of the reference trajectory in global coordinates, and that global information is not available in the particle data alone, so it would need to be specified in some other way.

The reference values of :math:`z=ct` and :math:`s` should be identical until reaching a bending element.
(If the lattice contains no bending elements, then they should coincide.)
Also, the reference value of :math:`ct` coincides with the value of :math:`s` in the ultrarelativistic limit.
More generally, the derivative :math:`ds/d(ct) = \beta`, where the relativistic :math:`\beta = \sqrt{1-\frac{1}{p_t^2}}`.


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
