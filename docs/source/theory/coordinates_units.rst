.. _theory-coordinates-and-units:

Coordinates and Units
=====================

Each particle in the beam is described at fixed :math:`s` by a set of 6 canonical phase space variables (x [m], px, y [m], py, t [m], pt).  Coordinates x and y denote the horizontal and
vertical displacement from the reference particle, respectively, and describe motion in the plane transverse to the velocity of the reference particle.  The longitudinal coordinate t
denotes the difference between the arrival time of the particle and the arrival time of the reference particle, multiplied by the speed of light :math:`c`.

The momenta conjugate to x, y, and t are denoted px, py, and pt, respectively.  These variables are normalized by the magnitude of the momentum of the reference particle, and are therefore dimensionless.
In a region of zero vector potential, for example, :math:`p_x = \Delta(\beta_x\gamma)/(\beta_0\gamma_0)`, where :math:`\beta_0` and :math:`\gamma_0` denote the relativistic
factors associated with the reference velocity.  In a region of zero scalar potential, pt denotes the deviation from the reference energy normalized by the design momentum
times the speed of light, so that :math:`p_t = -\Delta(\gamma)/(\beta_0\gamma_0)`.

Unlike particles within the beam, the reference particle is described by a set of 8 phase space variables (x [m], px, y [m], py, z [m], pz, t [m], pt) that are specified
in a global laboratory coordinate system (x,y,z).  The momenta of the reference particle are normalized by :math:`mc`, so that :math:`p_x=\beta_x\gamma`, etc.  A parameteric plot of
the reference trajectory variables (x,z) allows the user to view the global geometry of the accelerator structure (footprint).
