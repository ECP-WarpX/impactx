.. _theory-assumptions:

Assumptions
===========

This is an overview of physical assumptions implemented in the numerics of ImpactX.


Tracking and Lattice Optics
---------------------------

Tracking through lattice optics in ImpactX is performed by updating the canonical phase space variables (x,px,y,py,t,pt) using symplectic transport.
The elements supported currently fall into one of the following categories:

* **zero-length (thin) elements**, such as multipole kicks and coordinate transformations
* **ideal (thick) elements** using a hard-edge fringe field approximation, such as drifts, quadrupoles, and dipoles
* **soft-edge elements** described by :math:`s`-dependent, user-provided field data, such as RF cavities
* **ML surrogate models** using a trained neural network (not necessarily symplectic)

Transport may be performed using one of three possible levels of approximation to the underlying Hamiltonian:

* **linear transfer map (default):** obtained by expanding the Hamiltonian through terms of degree 2 in the deviation of phase space variables from those of the reference particle
* **chromatic or paraxial approximation:** obtained by expanding the Hamiltonian through terms of degree 2 in the transverse phase space variables, while retaining the nonlinear dependence on the energy variable pt
* **exact Hamiltonian:** obtained using the exact nonlinear Hamiltonian


Space Charge (Poisson Solver)
-----------------------------

  * **velocity spread:** when solving for space-charge effects, we assume that the relative spread of velocities of particles in the beam is negligible compared to the velocity of the reference particle, so that in the bunch frame (rest frame of the reference particle) particle velocities are nonrelativistic

* **electrostatic in the bunch frame:** we assume there are no retardation effects and we solve the Poisson equation in the bunch frame
