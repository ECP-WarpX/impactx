.. _theory-collective-beam-distribution-input:

Beam Distribution Input
=======================

.. figure:: phase_space_ellipse.svg
   :align: center
   :width: 75%
   :alt: phase space ellipse

   Phase space ellipse in the coordinate plane of position :math:`q` (realized as :math:`x`, :math:`y`, and :math:`t`) and associated conjugate momentum :math:`q` (realized as :math:`p_x`, :math:`p_y` and :math:`p_t`).

Particle beam user input in ImpactX can be done in two ways.

The first option is to characterize the distribution via the intersections of the phase space ellipse with the coordinate axes and the correlation terms of the canonical coordinate pairs.

.. math::

   \begin{align}
        \lambda_q &= \sqrt{\frac{\epsilon}{\gamma}} \\
        \lambda_p &= \sqrt{\frac{\epsilon}{\beta}} \\
        \mu_{qp} &= \frac{\alpha}{\sqrt{\beta \gamma}}
   \end{align}

The units are :math:`[\lambda_q] = \mathrm{m}`, :math:`[\lambda_p] = \mathrm{rad}`, and :math:`[\mu_{qp}] = 1`.
To convert between normalized and unnormalized emittance, use the relation :math:`\epsilon_\mathrm{n} = (\beta\gamma)_\mathrm{ref} \cdot \epsilon` which uses the momentum of the reference particle.
**Attention**: Here, :math:`(\beta\gamma)_\mathrm{ref}` are the Lorentz variables for the reference particle momentum and not the Courant-Snyder parameters.

The second option is to specify the distribution via the Courant-Snyder / Twiss parameters :math:`\alpha` and :math:`\beta` along with the unnormalized (geometric, 1-RMS) emittance :math:`\epsilon` for all the spatial coordinates.
Recall the Courant-Snyder relation :math:`\gamma\beta - \alpha^2 = 1` for conversion from :math:`\gamma` values to our input conventions.

Distribution Sampling and the Covariance Matrix
===============================================

In ImpactX, beam sampling is performed under the assumption that the initial beam distribution centroid (mean phase space vector) coincides with the phase space origin.  The covariance matrix :math:`\Sigma` is defined by :math:`\Sigma_{ij}=\langle{\zeta_i\zeta_j\rangle}`, where :math:`\zeta` denotes the vector of phase space coordinates, and indices :math:`i,j` specify the components of :math:`\zeta`.

Let :math:`P` denote a phase space probability density with unit covariance matrix (i.e., equal to the identity matrix).  To produce a phase space density with a target covariance matrix :math:`\Sigma`, we write :math:`\Sigma` in terms of its (lower) Cholesky decomposition as:

.. math::

   \begin{equation}
        \Sigma = LL^T,
   \end{equation}

where :math:`L` is a lower triangular matrix.

Define a beam distribution function :math:`f` by:

.. math::

   \begin{equation}
       f(\zeta)=\kappa P(L^{-1}\zeta),\quad\text{where}\quad \kappa=|\det L|^{-1}.
   \end{equation}

Then :math:`f` has the desired covariance matrix :math:`\Sigma`.  Samples from :math:`f` are obtained by sampling from :math:`P` and performing the linear transformation :math:`\zeta\mapsto L\zeta`.

Let :math:`P` above denote a 2D probability distribution that is radially symmetric, in the sense that:

.. math::

   \begin{equation}
        P(\zeta)=G(||\zeta||^2)=G(q^2+p^2),\quad\quad \zeta=(q,p)
   \end{equation}

Here :math:`q` denotes a position coordinate (e.g., :math:`x`, :math:`y`, or :math:`t`) and :math:`p` denotes the corresponding conjugate momentum.

Then the resulting distribution :math:`f` has 2D elliptical symmetry, in the sense that:

.. math::

   \begin{equation}
        f(\zeta)\propto P(L^{-1}\zeta)=G(||L^{-1}\zeta||^2)=G(\zeta^TS\zeta),\quad\quad S=\Sigma^{-1}.
   \end{equation}

The argument of :math:`G` is a quadratic form in :math:`(q,p)`, and it is convenient to express this quadratic form as:

.. math::

   \begin{equation}
        \zeta^TS\zeta = \frac{q^2}{\lambda_q^2} + 2\mu_{qp}\frac{qp}{\lambda_q\lambda_p}+\frac{p^2}{\lambda_p^2}=\frac{1}{\epsilon}\left(\gamma q^2+2\alpha qp + \beta p^2\right).
   \end{equation}

Here :math:`\alpha`, :math:`\beta`, and :math:`\gamma` denote the Courant-Snyder Twiss functions, and :math:`\epsilon` denotes the rms (unnormalized) emittance.

The associated covariance matrix may be written explicitly in terms of the above parameters as:

.. math::

   \begin{equation}
        \begin{pmatrix}
            \lambda_q & 0 \\
            0 & \lambda_p
        \end{pmatrix}
        \begin{pmatrix}
            1 & \mu_{qp} \\
            \mu_{qp} & 1
        \end{pmatrix}^{-1}
        \begin{pmatrix}
            \lambda_q & 0 \\
            0 & \lambda_p
        \end{pmatrix} = \epsilon
        \begin{pmatrix}
            \beta & -\alpha \\
           -\alpha & \gamma
        \end{pmatrix}.
   \end{equation}

Note:  In the special case that :math:`\mu_{qp}=0`, we have :math:`\lambda_q=\sigma_q` and :math:`\lambda_p=\sigma_p`, where :math:`\sigma_q=\langle{q^2\rangle}^{1/2}` and :math:`\sigma_p=\langle{p^2\rangle}^{1/2}`.
