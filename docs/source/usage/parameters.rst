.. _running-cpp-parameters:

Parameters: Inputs File
=======================

This documents how to use ImpactX with an input file (``impactx input_file.in``).

.. tip::

   If you enjoy AI/LLM/agentic workflows, see our :ref:`AI (LLM)-Assisted Input File Design <ai_input_design>` section, too.

.. tip::

   Input files use the AMReX `ParmParse <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`__ syntax.
   A `parser <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parser>`__) is used for the right-hand-side of all input parameters that consist of one or more integers or floats, so expressions like ``beam.kin_energy = "2.+1."``, ``beam.lambdaY = beam.lambdaX`` and/or using :ref:`user-defined constants <running-cpp-parameters-parser>` are accepted.


.. _running-cpp-parameters-mode:

Tracking Modes
--------------

.. pp:param:: algo.track
    :type: ``string``
    :default: ``particles``

    Mode that specifies how the beam is tracked:

    * ``particles`` (default): symplectic particle tracking
    * ``envelope``: beam envelope (covariance matrix) tracking, through linearized transport maps
    * ``reference_orbit``: only tracking of the reference particle orbit

    .. note::

       Our current ``envelope`` tracking implements ideal transfer maps, assuming always zero misalignments (translations).
       Element rotations are handled.
       Support for translations errors, non-zero envelope means, and feed-down effects in envelope tracking is in development.
       Until then, translations errors set on elements are silently ignored.

.. _running-cpp-parameters-particle:

Initial Beam Distributions
--------------------------

.. pp:param:: beam.npart
    :type: ``integer``

    Number of weighted simulation particles.

.. pp:param:: beam.units
    :type: ``string``

    Currently, only ``static`` is supported.

.. pp:param:: beam.kin_energy
    :type: ``float``
    :unit: MeV

    Beam kinetic energy.

.. pp:param:: beam.charge
    :type: ``float``
    :unit: C

    Bunch charge.

.. pp:param:: beam.current
    :type: ``float``
    :unit: A

    Beam current, used only if :pp:param:`algo.space_charge = 2D`.

.. pp:param:: beam.particle
    :type: ``string``

    Particle type: ``electron``, ``positron``, ``proton``, or ``Hminus``.
    For other species, please use the :ref:`Python API <usage-python>` or `open an issue <https://github.com/BLAST-ImpactX/impactx/issues>`__.

    .. dropdown:: Species Constants
       :color: light
       :icon: info
       :animate: fade-in-slide-down

       .. literalinclude:: ../../../src/particles/ReferenceParticle.H
          :language: cpp
          :dedent: 12
          :start-after: // [known_species]
          :end-before: // [/known_species]

.. pp:param:: beam.distribution
    :type: ``string``

    Indicates the initial distribution type.
    For additional information, consult the documentation on :ref:`theory-collective-beam-distribution-input`.
    For **all** except the ``thermal`` distribution we allow input in two forms (see the parameters below the list of distributions).

    The following distributions are available:

    * ``waterbag`` or ``waterbag_from_twiss`` for initial Waterbag distribution.

    * ``kurth6d`` or ``kurth6d_from_twiss`` for initial 6D Kurth distribution.

    * ``gaussian`` or ``gaussian_from_twiss`` for initial 6D Gaussian (normal) distribution.

      Optionally, the user may specify an independent cutoff in each phase plane (x,px), (y,py), and (t,pt).
      The cut is performed in normalized Courant-Snyder variables corresponding to the user-supplied second moments or Twiss functions.
      As a result, this is equivalent to a cut corresponding to the (linearized) action in each plane.

      Optional parameters:

      .. pp:param:: beam.cutX/Y/T
          :link_aliases: beam.cutX beam.cutY beam.cutT
          :type: ``float``
          :unit: dimensionless
          :default: ``0``

          Number of sigma at which to cut the distribution in (x,px) / (y,py) / (t,pt); ``0`` means no cut.

    * ``kvdist`` or ``kvdist_from_twiss`` for initial K-V distribution in the transverse plane.

      The distribution is uniform in t and Gaussian in pt.

    * ``kurth4d`` or ``kurth4d_from_twiss`` for initial 4D Kurth distribution in the transverse plane.

      The distribution is uniform in t and Gaussian in pt.

    * ``semigaussian`` or ``semigaussian_from_twiss`` for initial Semi-Gaussian distribution.

      The distribution is uniform within a cylinder in (x,y,z) and Gaussian in momenta (px,py,pt).

    * ``triangle`` or ``triangle_from_twiss`` a triangle distribution for laser-plasma acceleration related applications.

      A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).
      The transverse distribution is a 4D waterbag.

    * ``thermal`` for a 6D stationary thermal or bithermal distribution.

      This distribution type is described, for example in:
      R. D. Ryne et al., `"A Test Suite of Space-Charge Problems for Code Benchmarking" <https://accelconf.web.cern.ch/e04/PAPERS/WEPLT047.PDF>`__, in Proc. EPAC2004, Lucerne, Switzerland.
      C. E. Mitchell et al., `"ImpactX Modeling of Benchmark Tests for Space Charge Validation" <https://doi.org/10.18429/JACoW-HB2023-THBP44>`__, in Proc. HB2023, Geneva, Switzerland.

      Additional parameters:

      .. pp:param:: beam.k
          :type: ``float``
          :unit: 1/m

          External focusing strength.

      .. pp:param:: beam.kT
          :type: ``float``
          :unit: dimensionless

          Temperature of core population
          :math:`= < p_x^2 > = < p_y^2 >`, where all momenta are normalized by the reference momentum.

      .. pp:param:: beam.kT_halo
          :type: ``float``
          :unit: dimensionless
          :default: ``kT``

          Temperature of halo population.

      .. pp:param:: beam.normalize
          :type: ``float``
          :unit: dimensionless

          Normalizing constant for core population.

      .. pp:param:: beam.normalize_halo
          :type: ``float``
          :unit: dimensionless

          Normalizing constant for halo population.

      .. pp:param:: beam.halo
          :type: ``float``
          :unit: dimensionless

          Fraction of charge in halo.

For **all** except the ``thermal`` distribution, the phase space ellipse is specified in one of two forms.

#. Parameters that describe the **phase space ellipse and position-momentum correlations**:

   .. pp:param:: beam.lambdaX/Y/T
       :link_aliases: beam.lambdaX beam.lambdaY beam.lambdaT
       :type: ``float``
       :unit: m

       Phase space ellipse intersection with X / Y / T.
       The T value is normalized by multiplying with the speed of light *c*.

   .. pp:param:: beam.lambdaPx/Py/Pt
       :link_aliases: beam.lambdaPx beam.lambdaPy beam.lambdaPt
       :type: ``float``
       :unit: rad

       Phase space ellipse intersection with Px / Py / Pt.

   .. pp:param:: beam.muxpx/muypy/mutpt
       :link_aliases: beam.muxpx beam.muypy beam.mutpt
       :type: ``float``
       :unit: dimensionless
       :default: ``0``

       Correlation X-Px / Y-Py / T-Pt.

#. **Courant-Snyder (Twiss)** parameters.
   To enable input via Courant-Snyder (Twiss) parameters, add the suffix ``from_twiss`` to the name of the distribution.
   Use the following parameters to characterize it:

   .. pp:param:: beam.alphaX/Y/T
       :link_aliases: beam.alphaX beam.alphaY beam.alphaT
       :type: ``float``
       :unit: dimensionless
       :default: ``0``

       CS / Twiss :math:`\alpha` for X / Y / T.

   .. pp:param:: beam.betaX/Y/T
       :link_aliases: beam.betaX beam.betaY beam.betaT
       :type: ``float``
       :unit: m

       CS / Twiss :math:`\beta` for X / Y / T.

   .. pp:param:: beam.emittX/Y/T
       :link_aliases: beam.emittX beam.emittY beam.emittT
       :type: ``float``
       :unit: m*rad

       Geometric (unnormalized) emittance :math:`\epsilon` in X / Y / T.

Two additional (optional) sets of input parameters may be provided:

#. Parameters that describe the displacement of the **beam centroid** from the reference particle in phase space:

   .. pp:param:: beam.meanX/Y/T
       :link_aliases: beam.meanX beam.meanY beam.meanT
       :type: ``float``
       :unit: m

       Mean value of the X / Y / T coordinate.

   .. pp:param:: beam.meanPx/Py/Pt
       :link_aliases: beam.meanPx beam.meanPy beam.meanPt
       :type: ``float``
       :unit: rad

       Mean value of the Px / Py / Pt coordinate.

#. Parameters that describe correlations between (x,pt), (px,pt), (y,pt), and (py,pt), described by a nonzero **dispersion**:

   .. pp:param:: beam.dispX/Y
       :link_aliases: beam.dispX beam.dispY
       :type: ``float``
       :unit: m

       Beam-based horizontal / vertical dispersion.

   .. pp:param:: beam.dispPx/Py
       :link_aliases: beam.dispPx beam.dispPy
       :type: ``float``
       :unit: dimensionless

       Derivative of beam-based horizontal / vertical dispersion.

.. pp:param:: beam.bucket_length
    :type: ``float``
    :unit: m

    Length of the longitudinal particle domain (e.g., length of the RF bucket in z), optionally provided for the application of particle boundary conditions.

Initial Spin Distributions
--------------------------

The specification of an initial particle spin distribution is optional, and is required only if spin tracking is used.
The default distribution type is the von Mises-Fisher distribution, uniquely determined by the input polarization vector.
The polarization vector provided by the user must lie within the unit ball.

.. pp:param:: beam.polarization_x/y/z
    :link_aliases: beam.polarization_x beam.polarization_y beam.polarization_z
    :type: ``float``
    :unit: dimensionless

    Mean value of the spin vector x / y / z-component.

.. _running-cpp-parameters-lattice:

Lattice Elements
----------------

.. pp:param:: lattice.elements
    :type: ``list of strings``
    :optional:
    :default: no elements

    A list of names (one name per lattice element), in the order that they appear in the lattice.

.. pp:param:: lattice.periods
    :type: ``integer``
    :optional:
    :default: ``1``

    The number of periods to repeat the lattice.

.. pp:param:: lattice.reverse
    :type: ``boolean``
    :optional:
    :default: ``false``

    Reverse the list of elements in the lattice.
    If :pp:param:`lattice.reverse` and :pp:param:`lattice.periods` both appear, then ``reverse`` is applied before ``periods``.

.. pp:param:: lattice.nslice
    :type: ``integer``
    :optional:
    :default: ``1``

    A positive integer specifying the number of slices used for the application of
    space charge in all elements; overwritten by the per-element ``nslice`` parameter.

.. pp:param:: <element_name>.type
    :type: ``string``

    Indicates the element type for this lattice element. This should be one of the following.

Many elements take a transverse aperture via ``aperture_x`` and ``aperture_y``: the ``aperture``
collimator, and as a beam pipe most other elements.
They all follow the same convention. Each plane is bounded independently: a half-aperture of zero
or less removes the constraint in that plane only, while the other plane still cuts.
Bounding a single plane gives a jaw (slit) collimator, bounding both an iris; with only
``aperture_y`` set, a particle is lost when ``|y| > aperture_y`` at any ``x``, and the
``rectangular`` and ``elliptical`` shapes degenerate to the same slab.
The aperture is disabled entirely only if both planes are zero or less, which is the default.


``aperture``
^^^^^^^^^^^^

``aperture`` for a thin collimator element applying a transverse aperture boundary,
e.g. for an element ``<aperture_name>.type = aperture``.
This requires these additional parameters:

.. pp:param:: <aperture_name>.aperture_x/y
    :link_aliases: <aperture_name>.aperture_x <aperture_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical or rectangular).

    Note that with ``action = absorb`` a disabled plane makes the absorbing domain
    unbounded in that direction, so zeroing both planes absorbs the whole beam,
    while for ``action = transmit`` the element has no effect on the beam.

.. pp:param:: <aperture_name>.repeat_x/y
    :link_aliases: <aperture_name>.repeat_x <aperture_name>.repeat_y
    :type: ``float``
    :unit: m

    Horizontal / vertical period for repeated aperture masking (inactive by default).

.. pp:param:: <aperture_name>.shift_odd_x
    :type: ``boolean``

    For hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by ``repeat_x / 2``.
    Use alignment offsets ``dx``, ``dy`` to move the whole mask as needed.

.. pp:param:: <aperture_name>.shape
    :type: ``string``
    :default: ``rectangular``

    Shape of the aperture boundary: ``rectangular`` (default) or ``elliptical``.

.. pp:param:: <aperture_name>.action
    :type: ``string``
    :default: ``transmit``

    Action of the aperture domain: ``transmit`` (default) or ``absorb``.

.. pp:param:: <aperture_name>.dx/dy
    :link_aliases: <aperture_name>.dx <aperture_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <aperture_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.


``beam_monitor``
^^^^^^^^^^^^^^^^

``beam_monitor`` a beam monitor, writing all beam particles at fixed ``s`` to openPMD files,
e.g. for an element ``<monitor_name>.type = beam_monitor``.
If the same element name is used multiple times, then an output series is created with multiple outputs.

.. pp:param:: <monitor_name>.name
    :type: ``string``
    :default: ``<monitor_name>``

    The output series name to use.
    By default, output is created under ``<diag.file_prefix>/openPMD/<monitor_name>.<backend>``.

.. pp:param:: <monitor_name>.backend
    :type: ``string``
    :default: ``default``

    `I/O backend <https://openpmd-api.readthedocs.io/en/latest/backends/overview.html>`_ for `openPMD <https://www.openPMD.org>`_ data dumps.
    ``bp4``/``bp5`` is the `ADIOS2 I/O library <https://csmd.ornl.gov/adios>`_, ``h5`` is the `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_, and ``json`` is a `simple text format <https://en.wikipedia.org/wiki/JSON>`_.
    ``json`` only works with serial/single-rank jobs.
    By default, the first available backend in the order given above is taken.

.. pp:param:: <monitor_name>.encoding
    :type: ``string``
    :default: ``g``

    openPMD `iteration encoding <https://openpmd-api.readthedocs.io/en/latest/usage/concepts.html#iteration-and-series>`__: (v)ariable based, (f)ile based, (g)roup based (default).
    Variable based is an experimental encoding and needs the ADIOS2 backend to store more than one iteration.
    openPMD-api 0.17 and newer read it with random access (``Access.read_only``) as well as with linear access (``Access.read_linear``); readers older than 0.17 see only a single iteration.

.. pp:param:: <monitor_name>.period_sample_intervals
    :type: ``integer``
    :default: ``1``

    For periodic lattice, only output every Nth period (turn).
    By default, diagnostics are returned every cycle.

.. pp:param:: <monitor_name>.particles
    :type: ``boolean``
    :default: ``true``

    Output the individual particles of the beam.
    When set to ``false``, the beam monitor becomes a fast, small diagnostic that stores only per-output attributes (reference particle and, if enabled, beam moments).

.. pp:param:: <monitor_name>.beam_moments
    :type: ``boolean``
    :default: ``true``

    Compute and output the beam moments (reduced beam characteristics) of the whole bunch as openPMD attributes on the ``beam`` particle species.

.. pp:param:: <monitor_name>.nonlinear_lens_invariants
    :type: ``boolean``
    :default: ``false``

    Compute and output the invariants H and I within the nonlinear magnetic insert element (see: ``nonlinear_lens``).
    Invariants associated with the nonlinear magnetic insert described by V. Danilov and S. Nagaitsev, PRSTAB 13, 084002 (2010), Sect. V.A.
    Requires ``<monitor_name>.particles = true``.

    .. pp:param:: <monitor_name>.alpha
        :type: ``float``
        :unit: dimensionless

        Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.

    .. pp:param:: <monitor_name>.beta
        :type: ``float``
        :unit: m

        Twiss beta of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.

    .. pp:param:: <monitor_name>.tn
        :type: ``float``
        :unit: dimensionless

        Dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I.

    .. pp:param:: <monitor_name>.cn
        :type: ``float``
        :unit: m^(1/2)

        Scale factor of the IOTA nonlinear magnetic insert element used for computing H and I.


``buncher``
^^^^^^^^^^^

``buncher`` for a short RF cavity (linear) bunching element,
e.g. for an element ``<buncher_name>.type = buncher``.
This requires these additional parameters:

.. pp:param:: <buncher_name>.V
    :type: ``float``
    :unit: dimensionless

    Normalized voltage drop across the cavity
    = (maximum voltage drop in Volts) / (speed of light in m/s * magnetic rigidity in T-m).

.. pp:param:: <buncher_name>.k
    :type: ``float``
    :unit: 1/m

    The RF wavenumber = 2*pi/(RF wavelength in m).

.. pp:param:: <buncher_name>.dx/dy
    :link_aliases: <buncher_name>.dx <buncher_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <buncher_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <buncher_name>.aperture_x/y
    :link_aliases: <buncher_name>.aperture_x <buncher_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``cfbend``
^^^^^^^^^^

``cfbend`` for a combined function bending magnet,
e.g. for an element ``<cfbend_name>.type = cfbend``.
This requires these additional parameters:

.. pp:param:: <cfbend_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <cfbend_name>.rc
    :type: ``float``
    :unit: m

    The bend radius.

.. pp:param:: <cfbend_name>.k
    :type: ``float``
    :unit: 1/m^2

    The quadrupole strength = (magnetic field gradient in T/m) / (magnetic rigidity in T-m).

    * ``k > 0`` horizontal focusing
    * ``k < 0`` horizontal defocusing

.. pp:param:: <cfbend_name>.dx/dy
    :link_aliases: <cfbend_name>.dx <cfbend_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <cfbend_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <cfbend_name>.aperture_x/y
    :link_aliases: <cfbend_name>.aperture_x <cfbend_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <cfbend_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``cfbend_exact``
^^^^^^^^^^^^^^^^

A thick combined-function dipole magnet using the exact relativistic Hamiltonian, including all kinematic nonlinearities.
The user must provide arrays containing normal and skew multipole coefficients, which can be specified up to decapole order.
The multipole coefficients are defined in the curvilinear coordinate system defined by the nominal reference trajectory.
For definitions of the coordinate system and (curvilinear) multipole coefficients we follow:

T. Zolkin, Phys. Rev. Accel. Beams 20, 043501 (2017), `DOI:10.1103/PhysRevAccelBeams.20.043501 <https://link.aps.org/doi/10.1103/PhysRevAccelBeams.20.043501>`__

The coefficients must appear in the following sequence:

dipole, quadrupole, sextupole, octupole, etc...

Particle tracking is performed using symplectic integration based on the Hamiltonian splitting :math:`H = H_1 + H_2`.
Here :math:`H_1` is the exact nonlinear Hamiltonian for a sector bend (including the kinematic square root),
and :math:`H_2` is the term containing the vector potential, which is a superposition of multipole contributions.

The vector potential is obtained from Table XI of the above-cited reference.

This element is defined via ``<cfbend_exact_name>.type = cfbend_exact`` and requires these additional parameters:

.. pp:param:: <cfbend_exact_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <cfbend_exact_name>.k_normal
    :type: ``array of float``

    Array of normal multipole coefficients (in meters^(-m) OR in T/meters^(m-1) for :math:`m=1,2,3,...`).

.. pp:param:: <cfbend_exact_name>.k_skew
    :type: ``array of float``

    Array of skew multipole coefficients (in meters^(-m) OR in T/meters^(m-1) for :math:`m=1,2,3,...`).

.. pp:param:: <cfbend_exact_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units for the multipole coefficients.
    By default, the multipole coefficients are normalized by magnetic rigidity. Use ``unit=1`` to specify using SI units.

.. pp:param:: <cfbend_exact_name>.dx/dy
    :link_aliases: <cfbend_exact_name>.dx <cfbend_exact_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <cfbend_exact_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <cfbend_exact_name>.aperture_x/y
    :link_aliases: <cfbend_exact_name>.aperture_x <cfbend_exact_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <cfbend_exact_name>.int_order
    :type: ``integer``
    :default: ``2``

    The order used for symplectic integration (2, 4, or 6).

.. pp:param:: <cfbend_exact_name>.mapsteps
    :type: ``integer``
    :default: ``10``

    Number of integration steps per slice used for symplectic integration.

.. pp:param:: <cfbend_exact_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``constf``
^^^^^^^^^^

``constf`` for a constant focusing element,
e.g. for an element ``<constf_name>.type = constf``.
This requires these additional parameters:

.. pp:param:: <constf_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <constf_name>.kx
    :type: ``float``
    :unit: 1/m

    The horizontal focusing strength.

.. pp:param:: <constf_name>.ky
    :type: ``float``
    :unit: 1/m

    The vertical focusing strength.

.. pp:param:: <constf_name>.kt
    :type: ``float``
    :unit: 1/m

    The longitudinal focusing strength.

.. pp:param:: <constf_name>.dx/dy
    :link_aliases: <constf_name>.dx <constf_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <constf_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <constf_name>.aperture_x/y
    :link_aliases: <constf_name>.aperture_x <constf_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <constf_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``dipedge``
^^^^^^^^^^^

``dipedge`` for dipole edge focusing. The model here is based on:

K. Hwang and S. Y. Lee, "Dipole fringe field map for compact synchrotrons," Phys. Rev. Accel. Beams 18, 122401 (2015)

as represented in the explicit, symplectic form provided in:

C. Mitchell and K. Hwang, "Explicit symplectic representations of nonlinear dipole fringe field maps," in Proc. NAPAC2025, TUP040, Sacramento, CA, 2025

Here, ``g`` denotes the magnetic gap, which is a length scale that sets the rate of decay of the fringe field.  The values ``K0`` - ``K6`` denote dimensionless field integrals, describing the shape of the fringe field, as defined in eqs. (28-34) of the first reference above.  In particular, ``K2`` is the well-known fringe field parameter denoted ``FINT`` in MAD-X.  The default values of the field integrals ``K0`` - ``K6`` are those given in eq. (52), corresponding to a ``tanh`` (i.e. logistic) field profile.

If ``model = "linear"``, then the linearized map is used.  This model is identical to:

* K. L. Brown, SLAC Report No. 75 (1982)

when expanded to first order in ``g/rc`` (gap / radius of curvature).

By comparison, note that the MAD-X DIPEDGE element uses as input the half-gap ``HGAP = g/2``, and sets the default value ``FINT = 0`` (while the corresponding default value of ``K2`` is set to 1).

Note that the nonlinear model includes a nonzero horizontal translation (depending on the field integral values) that is present even for a particle that begins on the ideal "hard-edge" reference trajectory.

For a beam, this will result in a centroid offset that will produce centroid oscillations in the  downstream beamline.
In practice, this can be avoided by aligning the downstream elements with the true horizontal position (after including the effect of the fringe field).
To model this correction, we allow two options in the dipedge model:

* the option ``modify_ref_part = False`` (default), in which the shift due to the fringe field is applied to each beam particle phase space vector but not to the reference particle phase space vector --
  this model makes sense if the shift due to the fringe field is not considered in the baseline design, so that downstream elements are aligned with the "idealized" reference trajectory

* the option ``modify_ref_part = True`` in which the shift due to the fringe field is applied to the reference particle phase space vector, but not to the beam particle phase space vector --
  this model makes sense if the shift due to the fringe field is considered as part of the baseline design, so that downstream elements are aligned with the "shifted" reference trajectory

This element is defined via ``<dipedge_name>.type = dipedge`` and requires these additional parameters:

.. pp:param:: <dipedge_name>.psi
    :type: ``float``
    :unit: rad

    The pole face rotation angle.

.. pp:param:: <dipedge_name>.rc
    :type: ``float``
    :unit: m

    The bend radius.

.. pp:param:: <dipedge_name>.g
    :type: ``float``
    :unit: m

    The full magnetic gap size.

.. pp:param:: <dipedge_name>.R
    :type: ``float``
    :unit: m
    :default: ``1``

    Scale length for the field integrals.

.. pp:param:: <dipedge_name>.K0
    :type: ``float``
    :unit: dimensionless
    :default: ``pi**2/6``

    Normalized field integral for fringe field.

.. pp:param:: <dipedge_name>.K1
    :type: ``float``
    :unit: dimensionless
    :default: ``0``

    Normalized field integral for fringe field.

.. pp:param:: <dipedge_name>.K2
    :type: ``float``
    :unit: dimensionless
    :default: ``1``

    Normalized field integral for fringe field (FINT).

.. pp:param:: <dipedge_name>.K3
    :type: ``float``
    :unit: dimensionless
    :default: ``1/6``

    Normalized field integral for fringe field.

.. pp:param:: <dipedge_name>.K4
    :type: ``float``
    :unit: dimensionless
    :default: ``0``

    Normalized field integral for fringe field.

.. pp:param:: <dipedge_name>.K5
    :type: ``float``
    :unit: dimensionless
    :default: ``0``

    Normalized field integral for fringe field.

.. pp:param:: <dipedge_name>.K6
    :type: ``float``
    :unit: dimensionless
    :default: ``0``

    Normalized field integral for fringe field.

.. pp:param:: <dipedge_name>.model
    :type: ``string``
    :default: ``linear``

    The fringe field model: ``linear`` (default) or ``nonlinear``.

.. pp:param:: <dipedge_name>.location
    :type: ``string``
    :default: ``entry``

    The fringe field edge location: ``entry`` (default) or ``exit``.

.. pp:param:: <dipedge_name>.modify_ref_part
    :type: ``boolean``
    :default: ``false``

    Apply fringe field to the reference particle, ``true`` or ``false`` (default).

.. pp:param:: <dipedge_name>.dx/dy
    :link_aliases: <dipedge_name>.dx <dipedge_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <dipedge_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <dipedge_name>.aperture_x/y
    :link_aliases: <dipedge_name>.aperture_x <dipedge_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``drift``
^^^^^^^^^

``drift`` for a free drift,
e.g. for an element ``<drift_name>.type = drift``.
This requires these additional parameters:

.. pp:param:: <drift_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <drift_name>.dx/dy
    :link_aliases: <drift_name>.dx <drift_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <drift_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <drift_name>.aperture_x/y
    :link_aliases: <drift_name>.aperture_x <drift_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <drift_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``drift_chromatic``
^^^^^^^^^^^^^^^^^^^

``drift_chromatic`` for a free drift, with chromatic effects included,
e.g. for an element ``<drift_chromatic_name>.type = drift_chromatic``.
The Hamiltonian is expanded through second order in the transverse variables (x,px,y,py), with the exact pt dependence retained.
This requires these additional parameters:

.. pp:param:: <drift_chromatic_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <drift_chromatic_name>.dx/dy
    :link_aliases: <drift_chromatic_name>.dx <drift_chromatic_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <drift_chromatic_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <drift_chromatic_name>.aperture_x/y
    :link_aliases: <drift_chromatic_name>.aperture_x <drift_chromatic_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <drift_chromatic_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``drift_exact``
^^^^^^^^^^^^^^^

``drift_exact`` for a free drift, using the exact nonlinear map,
e.g. for an element ``<drift_exact_name>.type = drift_exact``.
This requires these additional parameters:

.. pp:param:: <drift_exact_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <drift_exact_name>.dx/dy
    :link_aliases: <drift_exact_name>.dx <drift_exact_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <drift_exact_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <drift_exact_name>.aperture_x/y
    :link_aliases: <drift_exact_name>.aperture_x <drift_exact_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <drift_exact_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``kicker``
^^^^^^^^^^

``kicker`` for a thin transverse kicker,
e.g. for an element ``<kicker_name>.type = kicker``.
This requires these additional parameters:

.. pp:param:: <kicker_name>.xkick/ykick
    :link_aliases: <kicker_name>.xkick <kicker_name>.ykick
    :type: ``float``
    :unit: dimensionless OR T-m

    The horizontal / vertical kick strength.

.. pp:param:: <kicker_name>.unit
    :type: ``string``
    :default: ``dimensionless``

    Specification of units: ``dimensionless`` (default, in units of the magnetic rigidity of the reference particle) or ``T-m``.

.. pp:param:: <kicker_name>.dx/dy
    :link_aliases: <kicker_name>.dx <kicker_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <kicker_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <kicker_name>.aperture_x/y
    :link_aliases: <kicker_name>.aperture_x <kicker_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``line``
^^^^^^^^

``line`` a sub-lattice (line) of elements to append to the lattice,
e.g. for an element ``<line_name>.type = line``.

.. pp:param:: <line_name>.elements
    :type: ``list of strings``
    :optional:
    :default: no elements

    A list of names (one name per lattice element), in the order that they appear in the lattice.

.. pp:param:: <line_name>.reverse
    :type: ``boolean``
    :optional:
    :default: ``false``

    Reverse the list of elements in the line before appending to the lattice.

.. pp:param:: <line_name>.repeat
    :type: ``integer``
    :optional:
    :default: ``1``

    Repeat the line multiple times before appending to the lattice.
    Note: If :pp:param:`<line_name>.reverse` and :pp:param:`<line_name>.repeat` both appear, then ``reverse`` is applied before ``repeat``.


``linear_map``
^^^^^^^^^^^^^^

``linear_map`` for a custom, linear transport matrix.

The matrix elements :math:`R(i,j)` are indexed beginning with 1, so that :math:`i,j=1,2,3,4,5,6`.
The transport matrix :math:`R` is defaulted to the identity matrix, so only matrix entries that differ from that need to be specified.

The matrix :math:`R` multiplies the phase space vector :math:`(x,px,y,py,t,pt)`, where coordinates :math:`(x,y,t)` have units of m
and momenta :math:`(px,py,pt)` are dimensionless.  So, for example, :math:`R(1,1)` is dimensionless, and :math:`R(1,2)` has units of m.

The internal tracking methods used by ImpactX are symplectic.  However, if a user-defined linear map :math:`R` is provided, it is up to the user to ensure that the matrix :math:`R` is symplectic.  Otherwise, this condition may be violated.

.. note::

   This element cannot be sliced, so a collective effect kick (space charge, CSR, ISR) is applied as ``K(ds/2) M(ds) K(ds/2)``: half a kick, the element map, then the other half, at two collective solves per element (see :pp:param:`algo.strang_split`).
   Model it as several shorter elements to resolve collective effects along its length.

This element is defined via ``<linear_map_name>.type = linear_map`` and requires these additional parameters:

.. pp:param:: <linear_map_name>.R(i,j)
    :type: ``float``

    Matrix entries: a 1-indexed, 6x6, linear transport map to multiply with the phase space vector :math:`(x,p_x,y,p_y,t,p_t)`.

.. pp:param:: <linear_map_name>.ds
    :type: ``float``
    :unit: m
    :default: ``0``

    Length associated with a user-defined linear element.

.. pp:param:: <linear_map_name>.dx/dy
    :link_aliases: <linear_map_name>.dx <linear_map_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <linear_map_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <linear_map_name>.aperture_x/y
    :link_aliases: <linear_map_name>.aperture_x <linear_map_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``multipole``
^^^^^^^^^^^^^

``multipole`` for a thin multipole element,
e.g. for an element ``<multipole_name>.type = multipole``.
This requires these additional parameters:

.. pp:param:: <multipole_name>.multipole
    :type: ``integer``
    :unit: dimensionless

    Order of multipole: (m = 1) dipole, (m = 2) quadrupole, (m = 3) sextupole, etc.

.. pp:param:: <multipole_name>.k_normal
    :type: ``float``

    Integrated normal multipole coefficient (MAD-X convention), in meters^(-m+1)
    = ds * 1/(magnetic rigidity in T-m) * (derivative of order :math:`m-1` of :math:`B_y` with respect to :math:`x`).

.. pp:param:: <multipole_name>.k_skew
    :type: ``float``

    Integrated skew multipole strength (MAD-X convention), in 1/meters^(-m+1).

.. pp:param:: <multipole_name>.dx/dy
    :link_aliases: <multipole_name>.dx <multipole_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <multipole_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <multipole_name>.aperture_x/y
    :link_aliases: <multipole_name>.aperture_x <multipole_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``multipole_exact``
^^^^^^^^^^^^^^^^^^^

``multipole_exact`` for a thick multipole magnet using the exact relativistic Hamiltonian, including all kinematic nonlinearities.
The user must provide arrays containing normal and skew multipole coefficients, which can be specified up to arbitrarily high order.
The fields are assumed to be uniform along the longitudinal beamline coordinate (hard-edge model).
The coefficients must appear in the following sequence:

dipole, quadrupole, sextupole, octupole, etc...

(Note: Dipole coefficients are currently ignored, and will be supported in a separate combined-function dipole element.)

Particle tracking is performed using symplectic integration based on the Hamiltonian splitting :math:`H = H_1 + H_2`.
Here :math:`H_1` is the nonlinear Hamiltonian for a drift (including the kinematic square root),
and :math:`H_2` is the term containing the vector potential, which is a superposition of multipole contributions.

This element is defined via ``<multipole_exact_name>.type = multipole_exact`` and requires these additional parameters:

.. pp:param:: <multipole_exact_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <multipole_exact_name>.k_normal
    :type: ``array of float``

    Array of normal multipole coefficients (in meters^(-m) OR in T/meters^(m-1) for :math:`m=1,2,3,...`).

.. pp:param:: <multipole_exact_name>.k_skew
    :type: ``array of float``

    Array of skew multipole coefficients (in meters^(-m) OR in T/meters^(m-1) for :math:`m=1,2,3,...`).

.. pp:param:: <multipole_exact_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units for the multipole coefficients.
    By default, the multipole coefficients are normalized by magnetic rigidity. Use ``unit=1`` to specify using SI units.

.. pp:param:: <multipole_exact_name>.dx/dy
    :link_aliases: <multipole_exact_name>.dx <multipole_exact_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <multipole_exact_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <multipole_exact_name>.aperture_x/y
    :link_aliases: <multipole_exact_name>.aperture_x <multipole_exact_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <multipole_exact_name>.int_order
    :type: ``integer``
    :default: ``2``

    The order used for symplectic integration (2, 4, or 6).

.. pp:param:: <multipole_exact_name>.mapsteps
    :type: ``integer``
    :default: ``10``

    Number of integration steps per slice used for symplectic integration.

.. pp:param:: <multipole_exact_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``nonlinear_lens``
^^^^^^^^^^^^^^^^^^

``nonlinear_lens`` for a thin IOTA nonlinear lens element,
e.g. for an element ``<nonlinear_lens_name>.type = nonlinear_lens``.
This requires these additional parameters:

.. pp:param:: <nonlinear_lens_name>.knll
    :type: ``float``
    :unit: m

    Integrated strength of the lens segment (MAD-X convention)
    = dimensionless lens strength * c parameter**2 * length / Twiss beta.

.. pp:param:: <nonlinear_lens_name>.cnll
    :type: ``float``
    :unit: m

    Distance of the singularities from the origin (MAD-X convention)
    = c parameter * sqrt(Twiss beta).

.. pp:param:: <nonlinear_lens_name>.dx/dy
    :link_aliases: <nonlinear_lens_name>.dx <nonlinear_lens_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <nonlinear_lens_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <nonlinear_lens_name>.aperture_x/y
    :link_aliases: <nonlinear_lens_name>.aperture_x <nonlinear_lens_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``plane_xyrotation``
^^^^^^^^^^^^^^^^^^^^

``plane_xyrotation`` for a rotation in the x-y plane (i.e., about the reference velocity vector),
e.g. for an element ``<plane_xyrotation_name>.type = plane_xyrotation``.
This requires these additional parameters:

.. pp:param:: <plane_xyrotation_name>.angle
    :type: ``float``
    :unit: degree

    Nominal angle of rotation.

.. pp:param:: <plane_xyrotation_name>.dx/dy
    :link_aliases: <plane_xyrotation_name>.dx <plane_xyrotation_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <plane_xyrotation_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.


``plasma_lens_chromatic``
^^^^^^^^^^^^^^^^^^^^^^^^^

``plasma_lens_chromatic`` for an active cylindrically-symmetric plasma lens, with chromatic effects included,
e.g. for an element ``<plasma_lens_chromatic_name>.type = plasma_lens_chromatic``.
The Hamiltonian is expanded through second order in the transverse variables :math:`(x,p_x,y,p_y)`, with the exact :math:`p_t` dependence retained.
This requires these additional parameters:

.. pp:param:: <plasma_lens_chromatic_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <plasma_lens_chromatic_name>.k
    :type: ``float``
    :unit: 1/m^2 OR T/m

    The plasma lens focusing strength
    = (azimuthal magnetic field gradient in T/m) / (magnetic rigidity in T-m) - if ``unit = 0``,

    OR = azimuthal magnetic field gradient in T/m - if ``unit = 1``.

.. pp:param:: <plasma_lens_chromatic_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units.

.. pp:param:: <plasma_lens_chromatic_name>.dx/dy
    :link_aliases: <plasma_lens_chromatic_name>.dx <plasma_lens_chromatic_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <plasma_lens_chromatic_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <plasma_lens_chromatic_name>.aperture_x/y
    :link_aliases: <plasma_lens_chromatic_name>.aperture_x <plasma_lens_chromatic_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <plasma_lens_chromatic_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``polygon_aperture``
^^^^^^^^^^^^^^^^^^^^

``polygon_aperture`` for a thin collimator element applying a transverse polygon aperture boundary defined by :math:`(x,y)` coordinates
and optional radius below which all particles are transmitted, e.g. for an element ``<polygon_aperture_name>.type = polygon_aperture``.
The vertices must define a closed curve and be ordered in the counter-clockwise direction.
The first and last vertices must be identical. These parameters define the element:

.. pp:param:: <polygon_aperture_name>.vertices_x/y
    :link_aliases: <polygon_aperture_name>.vertices_x <polygon_aperture_name>.vertices_y
    :type: ``array of float``
    :unit: m

    Array of horizontal / vertical locations of aperture vertices.

.. pp:param:: <polygon_aperture_name>.min_radius2
    :type: ``float``
    :unit: m^2
    :default: ``0``

    Optional minimum radius-squared of a circle fully inscribed within the polygon. Particles with
    radius-squared less than this value are transmitted by the aperture and the polygon calculation is skipped.

.. pp:param:: <polygon_aperture_name>.repeat_x/y
    :link_aliases: <polygon_aperture_name>.repeat_x <polygon_aperture_name>.repeat_y
    :type: ``float``
    :unit: m

    Horizontal / vertical period for repeated aperture masking (inactive by default).

.. pp:param:: <polygon_aperture_name>.shift_odd_x
    :type: ``boolean``

    For hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by ``repeat_x / 2``.
    Use alignment offsets ``dx``, ``dy`` to move the whole mask as needed.

.. pp:param:: <polygon_aperture_name>.action
    :type: ``string``
    :default: ``transmit``

    Action of the aperture domain: ``transmit`` (default) or ``absorb``.

.. pp:param:: <polygon_aperture_name>.dx/dy
    :link_aliases: <polygon_aperture_name>.dx <polygon_aperture_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <polygon_aperture_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

``prot``
^^^^^^^^

``prot`` for an exact pole-face rotation in the x-z plane,
e.g. for an element ``<prot_name>.type = prot``.
This requires these additional parameters:

.. pp:param:: <prot_name>.phi_in
    :type: ``float``
    :unit: degree

    Angle of the reference particle with respect to the longitudinal (z) axis in the original frame.

.. pp:param:: <prot_name>.phi_out
    :type: ``float``
    :unit: degree

    Angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame.


``quad``
^^^^^^^^

``quad`` for a quadrupole,
e.g. for an element ``<quad_name>.type = quad``.
This requires these additional parameters:

.. pp:param:: <quad_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <quad_name>.k
    :type: ``float``
    :unit: 1/m^2

    The quadrupole strength = (magnetic field gradient in T/m) / (magnetic rigidity in T-m).

    * ``k > 0`` horizontal focusing
    * ``k < 0`` horizontal defocusing

.. pp:param:: <quad_name>.dx/dy
    :link_aliases: <quad_name>.dx <quad_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <quad_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <quad_name>.aperture_x/y
    :link_aliases: <quad_name>.aperture_x <quad_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <quad_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``quad_chromatic``
^^^^^^^^^^^^^^^^^^

``quad_chromatic`` for a Quadrupole magnet, with chromatic effects included,
e.g. for an element ``<quad_chromatic_name>.type = quad_chromatic``.
The Hamiltonian is expanded through second order in the transverse variables :math:`(x,p_x,y,p_y)`, with the exact :math:`p_t` dependence retained.
This requires these additional parameters:

.. pp:param:: <quad_chromatic_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <quad_chromatic_name>.k
    :type: ``float``
    :unit: 1/m^2 OR T/m

    The quadrupole strength = (magnetic field gradient in T/m) / (magnetic rigidity in T-m) - if ``unit = 0``,

    OR = magnetic field gradient in T/m - if ``unit = 1``.

    * ``k > 0`` horizontal focusing
    * ``k < 0`` horizontal defocusing

.. pp:param:: <quad_chromatic_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units.

.. pp:param:: <quad_chromatic_name>.dx/dy
    :link_aliases: <quad_chromatic_name>.dx <quad_chromatic_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <quad_chromatic_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <quad_chromatic_name>.aperture_x/y
    :link_aliases: <quad_chromatic_name>.aperture_x <quad_chromatic_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <quad_chromatic_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``quad_exact``
^^^^^^^^^^^^^^^^^^

``quad_exact`` for a Quadrupole magnet using the exact relativistic Hamiltonian, including all kinematic nonlinearities.
Particle tracking is performed using symplectic integration based on the Hamiltonian splitting :math:`H = H_1 + H_2`.
Here :math:`H_1` is the Hamiltonian for a linear quadrupole (containing all terms quadratic in the phase space variables),
and :math:`H_2` is the remainder (including the kinematic square root).  This suggested splitting appears for example in:

* D. L. Bruhwiler et al, in Proc. of EPAC 98, pp. 1171-1173 (1998).
* E. Forest, J. Phys. A: Math. Gen. 39, 5321 (2006).

This element is defined via ``<quad_exact_name>.type = quad_exact`` and requires these additional parameters:

.. pp:param:: <quad_exact_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <quad_exact_name>.k
    :type: ``float``
    :unit: 1/m^2 OR T/m

    The quadrupole strength = (magnetic field gradient in T/m) / (magnetic rigidity in T-m) - if ``unit = 0``,

    OR = magnetic field gradient in T/m - if ``unit = 1``.

    * ``k > 0`` horizontal focusing
    * ``k < 0`` horizontal defocusing

.. pp:param:: <quad_exact_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units.

.. pp:param:: <quad_exact_name>.dx/dy
    :link_aliases: <quad_exact_name>.dx <quad_exact_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <quad_exact_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <quad_exact_name>.aperture_x/y
    :link_aliases: <quad_exact_name>.aperture_x <quad_exact_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <quad_exact_name>.int_order
    :type: ``integer``
    :default: ``2``

    The order used for symplectic integration (2, 4, or 6).

.. pp:param:: <quad_exact_name>.mapsteps
    :type: ``integer``
    :default: ``10``

    Number of integration steps per slice used for symplectic integration.

.. pp:param:: <quad_exact_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``quadrupole_softedge``
^^^^^^^^^^^^^^^^^^^^^^^

``quadrupole_softedge`` for a soft-edge quadrupole.  See :ref:`Models of Soft-Edge Elements <theory-softedge-elements>`.

The units used for the on-axis quadrupole gradient are the same as those used for the quadrupole strength ``k`` in the element Quad.  For example, if the values used to
describe the on-axis profile (as specified in ``cos_coefficients``, ``sin_coefficients``) attain a peak on-axis value of 1, then the parameter
``gscale``, which multiplies this profile, specifies the peak value of the quadrupole field gradient on-axis, divided by the magnetic rigidity.  In this case, ``gscale``
has units of inverse meters squared.

This element is defined via ``<quadrupole_softedge_name>.type = quadrupole_softedge`` and requires these additional parameters:

.. pp:param:: <quadrupole_softedge_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <quadrupole_softedge_name>.gscale
    :type: ``float``
    :unit: 1/m^2

    Scaling factor for on-axis magnetic field gradient.

.. pp:param:: <quadrupole_softedge_name>.cos_coefficients
    :type: ``array of float``
    :optional:

    Cos coefficients in Fourier expansion of the on-axis field gradient.
    Default is a tanh fringe field model from `MaryLie 3.0 <http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__.

.. pp:param:: <quadrupole_softedge_name>.sin_coefficients
    :type: ``array of float``
    :optional:

    Sin coefficients in Fourier expansion of the on-axis field gradient.
    Default is a tanh fringe field model from `MaryLie 3.0 <http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__.

.. pp:param:: <quadrupole_softedge_name>.dx/dy
    :link_aliases: <quadrupole_softedge_name>.dx <quadrupole_softedge_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <quadrupole_softedge_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <quadrupole_softedge_name>.aperture_x/y
    :link_aliases: <quadrupole_softedge_name>.aperture_x <quadrupole_softedge_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <quadrupole_softedge_name>.mapsteps
    :type: ``integer``
    :default: ``10``

    Number of integration steps per slice used for map and reference particle push in applied fields.

.. pp:param:: <quadrupole_softedge_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.

``quadedge``
^^^^^^^^^^^^

``quadedge`` for quadrupole edge focusing.  This is a nonlinear symplectic map (derived from a third-order Lie generator), representing the effect
of quadrupole entry or exit fringe fields in the hard-edge limit. This is an explicit symplectification of the Lie map that appears in eq (28) of:
E. Forest and J. Milutinovic, Nucl. Instrum. and Methods in Phys. Res. A 269, 474-482 (1988).
This element is defined via ``<quadedge_name>.type = quadedge`` and requires these additional parameters:

.. pp:param:: <quadedge_name>.k
    :type: ``float``
    :unit: 1/m^2 OR T/m

    The quadrupole strength = (magnetic field gradient in T/m) / (magnetic rigidity in T-m) - if ``unit = 0``,

    OR = magnetic field gradient in T/m - if ``unit = 1``.

    * ``k > 0`` horizontal focusing
    * ``k < 0`` horizontal defocusing

.. pp:param:: <quadedge_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units.

.. pp:param:: <quadedge_name>.flag
    :type: ``string``
    :default: ``entry``

    Specification of edge location: ``entry`` (default) or ``exit``.

.. pp:param:: <quadedge_name>.dx/dy
    :link_aliases: <quadedge_name>.dx <quadedge_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <quadedge_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <quadedge_name>.aperture_x/y
    :link_aliases: <quadedge_name>.aperture_x <quadedge_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

``rfcavity``
^^^^^^^^^^^^

``rfcavity`` a radiofrequency cavity.  See :ref:`Models of Soft-Edge Elements <theory-softedge-elements>`.

The units used for the on-axis longitudinal electric field are described in the documentation of ``escale`` below.  For example, if the values used to
describe the on-axis electric field (as specified in ``cos_coefficients``, ``sin_coefficients``, or ``gradient_on_axis``) attain a peak on-axis value of 1, then the parameter
``escale``, which multiplies this profile, specifies the peak value of the longitudinal electric field gradient on-axis, divided by particle rest energy.  In this case,
``escale`` has units of inverse meters.

This element is defined via ``<rfcavity_name>.type = rfcavity`` and requires these additional parameters:

.. pp:param:: <rfcavity_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <rfcavity_name>.escale
    :type: ``float``
    :unit: 1/m

    Scaling factor for on-axis RF electric field
    = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV).

.. pp:param:: <rfcavity_name>.freq
    :type: ``float``
    :unit: Hz

    RF frequency.

.. pp:param:: <rfcavity_name>.phase
    :type: ``float``
    :unit: degree

    RF driven phase.

.. pp:param:: <rfcavity_name>.cos_coefficients
    :type: ``array of float``
    :optional:

    Cosine coefficients in Fourier expansion of on-axis electric field Ez.
    Default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__.

.. pp:param:: <rfcavity_name>.sin_coefficients
    :type: ``array of float``
    :optional:

    Sine coefficients in Fourier expansion of on-axis electric field Ez.
    Default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__.

.. pp:param:: <rfcavity_name>.dx/dy
    :link_aliases: <rfcavity_name>.dx <rfcavity_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <rfcavity_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <rfcavity_name>.aperture_x/y
    :link_aliases: <rfcavity_name>.aperture_x <rfcavity_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <rfcavity_name>.mapsteps
    :type: ``integer``
    :default: ``10``

    Number of integration steps per slice used for map and reference particle push in applied fields.

.. pp:param:: <rfcavity_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``sbend``
^^^^^^^^^

``sbend`` for a bending magnet,
e.g. for an element ``<sbend_name>.type = sbend``.
This requires these additional parameters:

.. pp:param:: <sbend_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <sbend_name>.rc
    :type: ``float``
    :unit: m

    The bend radius.

.. pp:param:: <sbend_name>.dx/dy
    :link_aliases: <sbend_name>.dx <sbend_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <sbend_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <sbend_name>.aperture_x/y
    :link_aliases: <sbend_name>.aperture_x <sbend_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <sbend_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``sbend_exact``
^^^^^^^^^^^^^^^

``sbend_exact`` for a bending magnet using the exact nonlinear map for the bend body. The map corresponds to the map described in:
D. L. Bruhwiler et al., in Proc. of EPAC 98, pp. 1171-1173 (1998), E. Forest et al., Part. Accel. 45, pp. 65-94 (1994).  The model
consists of a uniform bending field B_y with a hard edge.  Pole faces are normal to the entry and exit velocity of the reference
particle.
This element is defined via ``<sbend_exact_name>.type = sbend_exact`` and requires these additional parameters:

.. pp:param:: <sbend_exact_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <sbend_exact_name>.phi
    :type: ``float``
    :unit: degree

    The bend angle.

.. pp:param:: <sbend_exact_name>.B
    :type: ``float``
    :unit: T
    :default: ``0``

    The bend magnetic field; when ``B = 0`` (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B.

.. pp:param:: <sbend_exact_name>.dx/dy
    :link_aliases: <sbend_exact_name>.dx <sbend_exact_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <sbend_exact_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <sbend_exact_name>.aperture_x/y
    :link_aliases: <sbend_exact_name>.aperture_x <sbend_exact_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <sbend_exact_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``shortrf``
^^^^^^^^^^^

``shortrf`` for a short RF cavity element,
e.g. for an element ``<shortrf_name>.type = shortrf``.
This requires these additional parameters:

.. pp:param:: <shortrf_name>.V
    :type: ``float``
    :unit: dimensionless

    Normalized voltage drop across the cavity
    = (maximum energy gain in MeV) / (particle rest energy in MeV).

.. pp:param:: <shortrf_name>.freq
    :type: ``float``
    :unit: Hz

    The RF frequency.

.. pp:param:: <shortrf_name>.phase
    :type: ``float``
    :unit: degree

    The synchronous RF phase.

    * ``phase = 0``: maximum energy gain (on-crest)
    * ``phase = -90 deg``: zero energy gain for bunching
    * ``phase = 90 deg``: zero energy gain for debunching

.. pp:param:: <shortrf_name>.dx/dy
    :link_aliases: <shortrf_name>.dx <shortrf_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <shortrf_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <shortrf_name>.aperture_x/y
    :link_aliases: <shortrf_name>.aperture_x <shortrf_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``solenoid``
^^^^^^^^^^^^

``solenoid`` for an ideal hard-edge solenoid magnet,
e.g. for an element ``<solenoid_name>.type = solenoid``.
This requires these additional parameters:

.. pp:param:: <solenoid_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <solenoid_name>.ks
    :type: ``float``
    :unit: 1/m

    Solenoid strength in m^(-1) (MADX convention) = (magnetic field Bz in T) / (rigidity in T-m).

.. pp:param:: <solenoid_name>.dx/dy
    :link_aliases: <solenoid_name>.dx <solenoid_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <solenoid_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <solenoid_name>.aperture_x/y
    :link_aliases: <solenoid_name>.aperture_x <solenoid_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <solenoid_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``solenoid_softedge``
^^^^^^^^^^^^^^^^^^^^^

``solenoid_softedge`` for a soft-edge solenoid.  See :ref:`Models of Soft-Edge Elements <theory-softedge-elements>`.

The units used for the on-axis longitudinal magnetic field data are determined by the parameter ``unit``.  For example, if the values used to
describe the on-axis profile (as specified in ``cos_coefficients``, ``sin_coefficients``, or ``field_on_axis``) attain a peak on-axis value of 1, then the parameter
``bscale``, which multiplies this profile, specifies the peak value of the longitudinal magnetic field gradient on-axis.  If ``unit=0``, this is normalized by the magnetic
rigidity.

This element is defined via ``<solenoid_softedge_name>.type = solenoid_softedge`` and requires these additional parameters:

.. pp:param:: <solenoid_softedge_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <solenoid_softedge_name>.bscale
    :type: ``float``
    :unit: 1/m

    Scaling factor for on-axis longitudinal magnetic field
    = (magnetic field Bz in T) / (magnetic rigidity in T-m) - if ``unit = 0``,

    OR = magnetic field Bz in T - if ``unit = 1``.

.. pp:param:: <solenoid_softedge_name>.cos_coefficients
    :type: ``array of float``
    :optional:

    Cos coefficients in Fourier expansion of the on-axis magnetic field Bz.
    Default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__.

.. pp:param:: <solenoid_softedge_name>.sin_coefficients
    :type: ``array of float``
    :optional:

    Sin coefficients in Fourier expansion of the on-axis magnetic field Bz.
    Default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__.

.. pp:param:: <solenoid_softedge_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units.

.. pp:param:: <solenoid_softedge_name>.dx/dy
    :link_aliases: <solenoid_softedge_name>.dx <solenoid_softedge_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <solenoid_softedge_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <solenoid_softedge_name>.aperture_x/y
    :link_aliases: <solenoid_softedge_name>.aperture_x <solenoid_softedge_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <solenoid_softedge_name>.mapsteps
    :type: ``integer``
    :default: ``10``

    Number of integration steps per slice used for map and reference particle push in applied fields.

.. pp:param:: <solenoid_softedge_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


``source``
^^^^^^^^^^

``source`` for a particle source,
e.g. for an element ``<source_name>.type = source``.
Typically at the beginning of a beam line.

Currently, this only supports openPMD files from our ``beam_monitor``.

In `MPI-parallel <https://www.mpi-forum.org>`__ runs, reading the particles in the source file is distributed over all MPI ranks:
each of the :math:`N` ranks reads a different contiguous chunk of :math:`1/N`-th of the particles, independent of their position.
When ImpactX needs to sort particles spatially, it will redistribute them over MPI ranks automatically during tracking.

.. pp:param:: <source_name>.distribution
    :type: ``string``

    Distribution type of particles in the source. Currently, only ``openPMD`` is supported.

.. pp:param:: <source_name>.openpmd_path
    :type: ``string``

    Path to the openPMD series.

.. pp:param:: <source_name>.active_once
    :type: ``boolean``
    :default: ``true``

    Inject particles only for the first lattice period.

.. pp:param:: <source_name>.load_ref_particle
    :type: ``boolean``
    :default: ``true``

    Restore the reference particle from the species metadata of the openPMD file.

    The reference particle state (:math:`s`, positions, momenta, mass, charge and the gyromagnetic anomaly)
    is restored exactly from the iteration the particles are read from.
    In this case, the :ref:`beam.* input block <running-cpp-parameters-particle>` can be omitted entirely.
    Set to ``false`` to load all particles relative to a manually configured reference particle:
    the relative particle coordinates in the file (in particular the momenta :math:`p_x`, :math:`p_y`
    and :math:`p_t`, which are normalized by the reference particle momentum) are then interpreted
    with respect to the existing reference particle.
    The existing reference particle in that case needs to be manually configured before the tracking loop.
    This option acts during particle tracking, it has no effect in envelope or reference-particle-only tracking.

.. pp:param:: <source_name>.load_step
    :type: ``integer``
    :default: unset (if neither option is set, the last step in the file is loaded)

    Which step (iteration) to load from the openPMD series, selected by step number.

    This is the ImpactX step at which the ``beam_monitor`` wrote the beam, which is stored as the
    openPMD iteration in the file.
    These step numbers are usually not consecutive, because the global step counter also advances in
    the elements between two monitors: list them with, e.g., `openpmd-ls <https://openpmd-api.readthedocs.io/en/0.17.1/utilities/cli.html>`__ before selecting one.
    A step that is not in the file is an error, which lists the steps that are in it.
    A negative value is an error, too: the step numbers in a file are not negative, use
    ``load_step_index`` to count back from the last step.

    Set at most one of ``load_step`` and ``load_step_index``.

.. pp:param:: <source_name>.load_step_index
    :type: ``integer``
    :default: unset (if neither option is set, the last step in the file is loaded)

    Which step (iteration) to load from the openPMD series, selected by position in the file.

    ``0`` is the first step in the file and ``-1`` is the last, counting back from it as in Python,
    e.g. ``-2`` is the second to last step.
    Use this when the step numbers in the file are not known, e.g. ``-2`` to continue from the turn
    before the last one that a ring wrote.
    An index that reaches past either end of the file is an error, which lists the steps in it.

    Set at most one of ``load_step`` and ``load_step_index``.


``spin_map``
^^^^^^^^^^^^

``spin_map`` for a custom, user-specified spin map that acts on the spin 3-vector :math:`(s_x,s_y,s_z)`.

Spin maps are specified in the Lie-algebraic form:

.. math::

   \vec{s}_f = M(\zeta)\vec{s}_i,\quad\quad M(\zeta)=e^{v\cdot L}e^{A\Delta\zeta\cdot L}.

Here :math:`v` is a 3-vector that defines the axis and angle of rotation at the phase space design point, and :math:`A` is a 3x6 matrix that defines the spin-orbit coupling for particles not on the design orbit.
Also, :math:`\Delta\zeta=(x,p_x,y,p_y,t,p_t)` denotes the 6-vector of phase space variables as deviations from the design orbit. The quantities :math:`L_x`, :math:`L_y`, and :math:`L_z` are standard 3x3 matrices that define a basis for the Lie algebra :math:`so(3)`.

The vector components :math:`v(i)` and the matrix elements :math:`A(i,j)` are indexed beginning with 1, so that :math:`i=1,2,3` and :math:`j=1,2,3,4,5,6`.
The vector :math:`v` and the matrix :math:`A` are defaulted to zero, so only entries that differ from zero need to be specified.

The matrix :math:`A` multiplies the phase space vector :math:`(x,p_x,y,p_y,t,p_t)`, where coordinates :math:`(x,y,t)` have units of m
and momenta :math:`(p_x,p_y,p_t)` are dimensionless.  The three components output are dimensionless.  So, for example, :math:`A(1,1)` has units of 1/m, and :math:`A(1,2)` is dimensionless.
All three components of :math:`v` are dimensionless.

.. note::

   This element cannot be sliced, so a collective effect kick (space charge, CSR, ISR) is applied as ``K(ds/2) M(ds) K(ds/2)``: half a kick, the element map, then the other half, at two collective solves per element (see :pp:param:`algo.strang_split`).
   Model it as several shorter elements to resolve collective effects along its length.

This element is defined via ``<spin_map_name>.type = spin_map`` and requires these additional parameters:

.. pp:param:: <spin_map_name>.v(i)
    :type: ``float``

    Vector entries: a 1-indexed, 3x1, axis-angle vector that defines the spin rotation at the phase space design point.

.. pp:param:: <spin_map_name>.A(i,j)
    :type: ``float``

    Matrix entries: a 1-indexed, 3x6, spin-orbit coupling matrix to multiply with the phase space vector :math:`(x,p_x,y,p_y,t,p_t)` that defines the spin rotation for off-design particles.

.. pp:param:: <spin_map_name>.ds
    :type: ``float``
    :unit: m
    :default: ``0``

    Length associated with a user-defined spin map element.

.. pp:param:: <spin_map_name>.dx/dy
    :link_aliases: <spin_map_name>.dx <spin_map_name>.dy
    :type: ``float``
    :unit: m
    :default: ``0``

    Horizontal / vertical translation error.

.. pp:param:: <spin_map_name>.rotation
    :type: ``float``
    :unit: degree
    :default: ``0``

    Rotation error in the transverse plane.


``tapered_pl``
^^^^^^^^^^^^^^

``tapered_pl`` for a thin nonlinear plasma lens with transverse (horizontal) taper.

.. math::

   B_x = g \left( y + \frac{xy}{D_x} \right), \quad \quad B_y = -g \left(x + \frac{x^2 + y^2}{2 D_x} \right)

where :math:`g` is the (linear) field gradient in T/m and :math:`D_x` is the targeted horizontal dispersion in m.

This element is defined via ``<tapered_pl_name>.type = tapered_pl`` and requires these additional parameters:

.. pp:param:: <tapered_pl_name>.k
    :type: ``float``
    :unit: 1/m OR T

    The integrated plasma lens focusing strength
    = (length in m) * (magnetic field gradient :math:`g` in T/m) / (magnetic rigidity in T-m) - if ``unit = 0``,

    OR = (length in m) * (magnetic field gradient :math:`g` in T/m) - if ``unit = 1``.

.. pp:param:: <tapered_pl_name>.unit
    :type: ``integer``
    :default: ``0``

    Specification of units.

.. pp:param:: <tapered_pl_name>.taper
    :type: ``float``
    :unit: 1/m

    Horizontal taper parameter = 1 / (target horizontal dispersion :math:`D_x` in m).

.. pp:param:: <tapered_pl_name>.dx/dy
    :link_aliases: <tapered_pl_name>.dx <tapered_pl_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <tapered_pl_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <tapered_pl_name>.aperture_x/y
    :link_aliases: <tapered_pl_name>.aperture_x <tapered_pl_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``thin_dipole``
^^^^^^^^^^^^^^^

``thin_dipole`` for a thin dipole element,
e.g. for an element ``<thin_dipole_name>.type = thin_dipole``.
This requires these additional parameters:

.. pp:param:: <thin_dipole_name>.theta
    :type: ``float``
    :unit: degree

    Dipole bend angle.

.. pp:param:: <thin_dipole_name>.rc
    :type: ``float``
    :unit: m

    Effective radius of curvature.

.. pp:param:: <thin_dipole_name>.dx/dy
    :link_aliases: <thin_dipole_name>.dx <thin_dipole_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <thin_dipole_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <thin_dipole_name>.aperture_x/y
    :link_aliases: <thin_dipole_name>.aperture_x <thin_dipole_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).


``uniform_acc_chromatic``
^^^^^^^^^^^^^^^^^^^^^^^^^

``uniform_acc_chromatic`` for a region of uniform acceleration, with chromatic effects included,
e.g. for an element ``<uniform_acc_chromatic_name>.type = uniform_acc_chromatic``.
The Hamiltonian is expanded through second order in the transverse variables (x,px,y,py), with the exact pt dependence retained.
This requires these additional parameters:

.. pp:param:: <uniform_acc_chromatic_name>.ds
    :type: ``float``
    :unit: m

    The segment length.

.. pp:param:: <uniform_acc_chromatic_name>.ez
    :type: ``float``
    :unit: 1/m

    The electric field strength
    = (particle charge in C * electric field Ez in V/m) / (particle mass in kg * (speed of light in m/s)^2).

.. pp:param:: <uniform_acc_chromatic_name>.bz
    :type: ``float``
    :unit: 1/m

    The magnetic field strength
    = (particle charge in C * magnetic field Bz in T) / (particle mass in kg * speed of light in m/s).

.. pp:param:: <uniform_acc_chromatic_name>.dx/dy
    :link_aliases: <uniform_acc_chromatic_name>.dx <uniform_acc_chromatic_name>.dy
    :type: ``float``
    :unit: m

    Horizontal / vertical translation error.

.. pp:param:: <uniform_acc_chromatic_name>.rotation
    :type: ``float``
    :unit: degree

    Rotation error in the transverse plane.

.. pp:param:: <uniform_acc_chromatic_name>.aperture_x/y
    :link_aliases: <uniform_acc_chromatic_name>.aperture_x <uniform_acc_chromatic_name>.aperture_y
    :type: ``float``
    :unit: m

    Horizontal / vertical half-aperture (elliptical).

.. pp:param:: <uniform_acc_chromatic_name>.nslice
    :type: ``integer``
    :default: ``1``

    Number of slices used for the application of space charge.


.. _running-cpp-parameters-collective:

Collective Effects
------------------

.. pp:param:: algo.strang_split
    :type: ``boolean``
    :optional:
    :default: ``true``

    Whether to compose the collective effect kicks and the element transport in a second-order, time-symmetric Strang split.
    Per slice, this applies ``M(ds/2) K(ds) M(ds/2)``: half of the slice transport ``M``, the collective kick ``K`` at the slice midpoint, then the other half.
    The collective solve still runs once per slice, so the second order comes at the cost of the first-order composition and far fewer slices are needed for the same accuracy.

    Elements that cannot be subdivided (``Programmable``, ``LinearMap`` and ``SpinMap``) halve the kick instead, ``K(ds/2) M(ds) K(ds/2)``, which needs two collective solves per slice.

    Setting this to ``false`` composes kick and transport to first order, ``K(ds) M(ds)``, as most other beam dynamics codes do.
    This is useful to compare against those codes and to study the convergence with ``nslice``.

.. _running-cpp-parameters-collective-spacecharge:

Space Charge
^^^^^^^^^^^^

Space charge kicks are applied per slice of thick :ref:`lattice elements <running-cpp-parameters-lattice>`.
See there ``nslice`` option on lattice elements for slicing.

.. pp:param:: algo.space_charge
    :type: ``string``
    :optional:
    :default: ``false``

    The physical model of space charge used.

    ImpactX uses an AMReX grid of boxes to organize and parallelize space charge simulation domain.
    These boxes also contain a field mesh, if space charge calculations are enabled.

    Options:

    * ``false`` (default): space charge effects are not calculated.

    * ``2D``: Space charge forces are computed in the plane ``(x,y)`` transverse to the reference particle velocity, assuming the beam is long and unbunched.

    * ``2p5D``: Space charge forces are computed in the plane ``(x,y)`` transverse to the reference particle velocity, while the transverse space charge kicks are weighted by the
      longitudinal line density determined by charge deposition (2.5D model).  Longitudinal space charge kicks are determined by the derivative of the line charge density.

      This model is supported only in particle tracking mode (when :pp:param:`algo.track = particles`).

      This model supports the following sub-options:

      .. pp:param:: algo.space_charge.num_longitudinal_bins
          :type: ``integer``
          :default: ``100``

          The number of bins for longitudinal line density deposition.

      .. pp:param:: algo.space_charge.apply_longitudinal_kick
          :type: ``boolean``
          :default: ``true``

          Enable or disable the longitudinal space charge kick.

    * ``3D``: Space charge forces are computed in three dimensions, assuming the beam is bunched.

      When running in envelope mode (when :pp:param:`algo.track = envelope`), this model currently assumes that ``<xy> = <yt> = <tx> = 0``.

    * ``Gauss3D``: Calculate 3D space charge forces as if the beam was a Gaussian distribution.

    * ``Gauss2p5D``: Calculate 2.5D space charge forces as if the beam was a transverse Gaussian distribution.

      These models are supported only in particle tracking mode (when :pp:param:`algo.track = particles`).
      Ref.: J. Qiang, "Two-and-a-half dimensional symplectic space-charge solver", LBNL Report Number: LBNL-2001674 (2025).
      (This reference describes both 3D and 2.5D models.)

      This model supports the following sub-option:

      .. pp:param:: algo.space_charge.gauss_nint
          :type: ``integer``
          :default: ``101``

          Number of steps for computing the integrals (Eqs. 45-47 in the above paper).

      .. pp:param:: algo.space_charge.gauss_taylor_delta
          :type: ``float``
          :default: ``0.01``

          Initial integral region to avoid divergence of integrand at 0.

      .. pp:param:: algo.space_charge.gauss_long_scale
          :type: ``float``
          :default: in-situ :math:`1.103 \cdot \gamma \cdot \sigma_z`

          Longitudinal space charge scale for the Gauss2p5D space charge model.
          This is an approximation that only influences the longitudinal momentum (``pt``) kick.
          If not set, it defaults to :math:`1.103 \cdot \gamma \cdot \sigma_z`, estimated in-situ from the
          current reduced beam characteristics (with :math:`\sigma_z` the RMS bunch length), which is a
          typical value when comparing to a 3D model.  This value minimizes the L2-norm of the error in
          the on-axis field Ez in the case of a long 3D Gaussian bunch.
          In a perfectly conducting pipe with a long bunch, where the bunch length in
          the rest frame is significantly larger than the pipe radius, it is recommended
          to set the longitudinal scale length to the pipe radius in the input file.

      .. pp:param:: algo.space_charge.gauss_charge_z_bins
          :type: ``integer``
          :default: ``129``

          Number of bins for longitudinal line density deposition.

.. pp:param:: amr.n_cell
    :type: ``3 integers``
    :optional:
    :default: 1 `blocking_factor <https://amrex-codes.github.io/amrex/docs_html/GridCreation.html>`__ per MPI process

    The number of grid points along each direction (on the **coarsest level**).

    If this parameter is not set, the mesh is decomposed into exactly one box
    (grid) per MPI process. One box per process is the ideal decomposition for
    every simulation, including space charge calculations. The automatically
    chosen mesh size is only suitable for simulations that need no auxiliary
    field grids, e.g., without space charge calculations. Otherwise, set this
    parameter to control the field resolution.

    If it is set, the mesh is chopped into boxes of at
    most :pp:param:`amr.max_grid_size` cells per direction, independent of the
    number of processes.

.. pp:param:: amr.max_grid_size
    :type: ``integer, per MR level``
    :optional:
    :default: ``32`` on CPU, ``64`` on GPU

    The maximum number of grid points per box (`grid <https://amrex-codes.github.io/amrex/docs_html/GridCreation.html>`__) along each direction.
    The mesh on each mesh-refinement level is decomposed into boxes of at most this size.

    If :pp:param:`amr.n_cell` is not set, ImpactX overwrites this default with the `blocking factor <https://amrex-codes.github.io/amrex/docs_html/GridCreation.html>`__ (``8``), which creates exactly one box per MPI process.

    For best performance, aim for **exactly one box per parallel (MPI) process**, which is also one box per GPU on accelerated machines.
    More boxes per process cause significant overhead (particle redistribution, communication in field solves, extra kernel launches).
    Fewer boxes than processes leave processes idle.
    ImpactX issues a warning at the start of a run if the number of boxes does not match the number of processes.

.. pp:param:: amr.max_grid_size_x
    :type: ``integer, per MR level``
    :optional:
    :default: ``amr.max_grid_size``

    The maximum number of grid points per box along ``x``, overwriting :pp:param:`amr.max_grid_size`.
    Analogous parameters exist as ``amr.max_grid_size_y`` and ``amr.max_grid_size_z``.

    The per-direction parameters always win over the isotropic :pp:param:`amr.max_grid_size`, independent of the order in which the two are given.

.. pp:param:: amr.max_level
    :type: ``integer``
    :default: ``0``

    When using mesh refinement, the number of refinement levels that will be used.

    Use ``0`` in order to disable mesh refinement.

.. pp:param:: amr.ref_ratio
    :type: ``integer``
    :default: ``2``

    When using mesh refinement, this is the refinement ratio per level.
    With this option, all directions are fined by the same ratio.

.. pp:param:: amr.ref_ratio_vect
    :type: ``3 integers``

    Set per refined level (for x,y,z).
    When using mesh refinement, this can be used to set the refinement ratio per direction and level, relative to the previous level.

    Example: for three levels, a value of ``2 2 4 8 8 16`` refines the first level by 2-fold in x and y and 4-fold in z compared to the coarsest level (level 0/mother grid); compared to the first level, the second level is refined 8-fold in x and y and 16-fold in z.

.. note::

   Particles that move outside the simulation domain are removed.

.. pp:param:: geometry.dynamic_size
    :type: ``boolean``
    :optional:
    :default: ``true``

    Use dynamic (``true``) resizing of the field mesh, via :pp:param:`geometry.prob_relative`, or static sizing (``false``), via :pp:param:`geometry.prob_lo` / :pp:param:`geometry.prob_hi`.

.. pp:param:: geometry.prob_relative
    :type: ``array of float``
    :unit: dimensionless
    :optional:
    :default: ``3.0 1.0 1.0 ...``

    A positive ``float`` array with :pp:param:`amr.max_level` entries.
    By default, we dynamically extract the minimum and maximum of the particle positions in the beam.
    The field mesh spans, per direction, multiple times the maximum physical extent of beam particles, as given by this factor.
    The beam minimum and maximum extent are symmetrically padded by the mesh.
    For instance, ``1.2`` means the mesh will span 10% above and 10% below the beam;
    ``1.0`` means the beam is exactly covered with the mesh.

.. pp:param:: geometry.prob_lo/hi
    :link_aliases: geometry.prob_lo geometry.prob_hi
    :type: ``3 floats``
    :unit: m
    :optional:

    Required if :pp:param:`geometry.dynamic_size` is ``false``.
    The extent of the full simulation domain relative to the reference particle position.
    This can be used to explicitly size the simulation box and ignore :pp:param:`geometry.prob_relative`.

    This box is rectangular, and thus its extent is given here by the coordinates of the lower corner (:pp:param:`geometry.prob_lo`) and upper corner (:pp:param:`geometry.prob_hi`).
    The first axis of the coordinates is x and the last is z.

.. pp:param:: algo.particle_shape
    :type: ``integer``

    Allowed values: ``1``, ``2``, or ``3``.
    The order of the shape factors (splines) for the macro-particles along all spatial directions: ``1`` for linear, ``2`` for quadratic, ``3`` for cubic.
    Low-order shape factors result in faster simulations, but may lead to more noisy results.
    High-order shape factors are computationally more expensive, but may increase the overall accuracy of the results.
    For production runs it is generally safer to use high-order shape factors, such as cubic order.

.. pp:param:: algo.poisson_solver
    :type: ``string``
    :optional:
    :default: ``fft``

    The numerical solver to solve the Poisson equation when calculating space charge effects.
    Currently, the multigrid solver supports only 3D space charge.  The fft solver supports either 2D or 3D space charge.
    An additional `2.5D solver <https://github.com/BLAST-ImpactX/impactx/issues/401>`__ will be added in the near future.

    Options:

    * ``fft``: Poisson's equation is solved using an Integrated Green Function method (which requires FFT calculations).

      See these references for more details `Qiang et al. (2006) <https://doi.org/10.1103/PhysRevSTAB.9.044204>`__ (+ `Erratum <https://doi.org/10.1103/PhysRevSTAB.10.129901>`__).
      This requires the compilation flag ``-DImpactX_FFT=ON``.
      If mesh refinement (MR) is enabled, this FFT solver is used only on the coarsest level and a multi-grid solver is used on refined levels.
      The boundary conditions are assumed to be open.

    * ``multigrid``: Poisson's equation is solved using an iterative multigrid (MLMG) solver.

      See the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#>`__ for details of the MLMG solver.
      Field boundaries for MLMG space charge calculation are located at the outer ends of the field mesh.
      For the MLMG solver, we assume `Dirichlet boundary conditions <https://en.wikipedia.org/wiki/Dirichlet_boundary_condition>`__ with zero potential (a mirror charge).
      Thus, to emulate open boundaries, consider adding enough vacuum padding to the beam.

Multigrid-specific numerical options:

.. pp:param:: algo.mlmg_relative_tolerance
    :type: ``float``
    :optional:
    :default: ``1.e-7`` (DP) / ``1.e-4`` (SP)

    The relative precision with which the electrostatic space-charge fields should be calculated.
    More specifically, the space-charge fields are computed with an iterative Multi-Level Multi-Grid (MLMG) solver.
    This solver can fail to reach the default precision within a reasonable time.

.. pp:param:: algo.mlmg_absolute_tolerance
    :type: ``float``
    :optional:
    :default: ``0``

    A value of ``0`` means: ignored.
    The absolute tolerance with which the space-charge fields should be calculated in units of V/m^2.
    More specifically, the acceptable residual with which the solution can be considered converged.
    In general this should be left as the default, but in cases where the simulation state changes very
    little between steps it can occur that the initial guess for the MLMG solver is so close to the
    converged value that it fails to improve that solution sufficiently to reach the
    :pp:param:`algo.mlmg_relative_tolerance` value.

.. pp:param:: algo.mlmg_max_iters
    :type: ``integer``
    :optional:
    :default: ``100``

    Maximum number of iterations used for MLMG solver for space-charge fields calculation.
    In case if MLMG converges but fails to reach the desired self_fields_required_precision,
    this parameter may be increased.

.. pp:param:: algo.mlmg_verbosity
    :type: ``integer``
    :optional:
    :default: ``1``

    The verbosity used for MLMG solver for space-charge fields calculation.
    Currently MLMG solver looks for verbosity levels from 0-5.
    A higher number results in more verbose output.


.. _running-cpp-parameters-collective-csr:

Coherent Synchrotron Radiation (CSR)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CSR effects are included in the simulation for bend lattice elements such as ``Sbend`` and ``CFbend``.
These effects are critical in accurately modeling the wakefields generated due to the interaction of particles with the synchrotron radiation field generated by the beam during bending.
Currently, this is the 1D ultrarelativistic steady-state wakefield model (eq. 19 of
`E. L. Saldin et al., NIMA 398, p. 373-394 (1997), DOI:10.1016/S0168-9002(97)00822-X <https://doi.org/10.1016/S0168-9002(97)00822-X>`__).

.. pp:param:: algo.csr
    :type: ``boolean``
    :optional:
    :default: ``false``

    Whether to calculate CSR effects.
    CSR calculations involve several steps, including charge deposition, wakefield generation, and convolution, all of which are handled within the CSR bending process.

.. pp:param:: algo.csr_bins
    :type: ``integer``
    :optional:
    :default: ``150``

    The number of bins used for the CSR calculations along the longitudinal direction. Increasing the number of bins can lead to more accurate wakefield resolution at the cost of higher computational expense.

.. note::

   CSR effects are only calculated for lattice elements that include bending, such as ``Sbend``, ``ExactSbend`` and ``CFbend``.

   CSR effects require the compilation flag ``-DImpactX_FFT=ON``.

.. _running-cpp-parameters-collective-isr:


Incoherent Synchrotron Radiation (ISR)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ISR effects are included in the simulation for bend lattice elements such as ``Sbend`` and ``CFbend``, and are especially important for electron or positron bunches at high energy.
The effects of ISR include radiation reaction due to the stochastic emission of synchrotron radiation, resulting in mean energy loss and quantum excitation of the bunch.
The model is based on:

`F. Niel et al., Phys. Rev. E 97, 043209 (2018), DOI:10.1103/PhysRevE.97.043209 <https://doi.org/10.1103/PhysRevE.97.043209>`__

However, a Taylor expansion is used to evaluate the dependence on the quantum parameter :math:`\chi`.  When ``algo.isr_order = 1``, the model is equivalent to that described in:

`J. M. Jowett, "Introductory Statistical Mechanics for Electron Storage Rings", AIP Conf. Proc. 153, 864-970 (1987), DOI:10.1063/1.36374 <https://doi.org/10.1063/1.36374>`__

.. pp:param:: algo.isr
    :type: ``boolean``
    :optional:
    :default: ``false``

    Whether to calculate ISR effects.

.. pp:param:: algo.isr_order
    :type: ``integer``
    :optional:
    :default: ``1``

    The number of terms retained in the Taylor series for the functions :math:`g(\chi)` and :math:`h(\chi)` appearing in equations (25) and (41) describing quantum effects.

.. pp:param:: algo.isr_on_ref_part
    :type: ``boolean``
    :optional:
    :default: ``false``

    Flag specifying whether ISR is to be applied to the reference particle.  When :pp:param:`algo.isr_on_ref_part = false`, the reference particle does not lose energy due to radiation, and the
    mean energy of the beam particles will decrease.  This option is natural if the lattice optics, magnet settings, etc. are chosen without accounting for radiative energy loss.
    When :pp:param:`algo.isr_on_ref_part = true`, the reference particle does lose energy due to radiation, and little centroid evolution is expected in the beam particles.  This option is natural if the lattice optics, magnet settings, etc. are chosen to account for radiative energy loss.

.. note::

   ISR effects are only calculated for lattice elements that include bending, such as ``Sbend``, ``ExactSbend`` and ``CFbend``.


.. _running-cpp-parameters-particle-bc:

Particle Boundary Condition
---------------------------

The application of a non-trivial boundary condition for particles is currently supported only in the longitudinal direction (the local direction of motion as defined by the reference particle).
For systems involving bunches that are long relative to the local size of the RF bucket, it is often necessary to capture the effect of particle slippage between adjacent buckets.
To handle this effectively without tracking multiple bunches, a periodic particle boundary condition can be applied.

.. pp:param:: algo.particle_bc
    :type: ``string``
    :optional:
    :default: ``open``

    The type of particle boundary condition to be applied to the longitudinal coordinate, based on the value of parameter :pp:param:`beam.bucket_length`.  Options:

    * ``open`` (default): no action is applied at the boundary.

    * ``periodic``: each particle's longitudinal coordinate is treated modulo the value :pp:param:`beam.bucket_length` (phase wrapping).

    * ``absorbing``: a particle whose longitudinal coordinate falls outside the boundary is declared lost.

    * ``reflecting``: a particle whose longitudinal coordinate crosses the boundary is reflected about the boundary, with reversed longitudinal momentum.

      .. note::

         The implementation works through linear order in the phase space variables.
         If you have the need for a more precise implementation of reflecting boundaries, please `open an issue <https://github.com/BLAST-ImpactX/impactx/issues/new>`__.


.. _running-cpp-parameters-spin:

Spin Tracking
-------------

Spin tracking is performed by updating the particle spin 3-vector in the presence of each element's electromagnetic fields, using methods based on the Thomas-BMT equation.
By construction, all spin tracking methods rely on pushing particles using spin maps that lie in SO(3).  The algorithm implemented for each element is consistent with the algorithm
used for pushing the phase space vector.

Currently, the implementation of spin tracking is a work in progress, and this feature is not yet supported.

.. pp:param:: algo.spin
    :type: ``boolean``
    :optional:
    :default: ``false``

    Whether to track particle spin.

    Spin tracking uses the gyromagnetic anomaly of the reference particle, which is set together with the particle species (see :pp:param:`beam.particle`).


.. _running-cpp-parameters-parser:

Math parser and user-defined constants
--------------------------------------

The AMReX parser is used for the right-hand-side of all input parameters that consist of one or more integers or floats.
Thus, expressions like ``beam.alphaY = beam.alphaX`` and/or using user-defined constants or simple math operations are accepted.

Note that when multiple values are expected, the expressions are space delimited.
For integer input values, the expressions are evaluated as real numbers and the final result rounded to the nearest integer.
See `this section <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parser>`__ of the AMReX documentation for a complete list of functions supported by the math parser.


ImpactX constants
^^^^^^^^^^^^^^^^^

ImpactX will provide a few pre-defined constants, that can be used for any parameter that consists of one or more floats.

======== ===================
q_e      elementary charge
m_e      electron mass
m_p      proton mass
m_u      unified atomic mass unit (Dalton)
epsilon0 vacuum permittivity
mu0      vacuum permeability
clight   speed of light
pi       math constant pi
======== ===================


User-defined constants
^^^^^^^^^^^^^^^^^^^^^^

Users can define their own constants in the input file.
These constants can be used for any parameter that consists of one or more integers or floats.
User-defined constant names can contain only letters, numbers and the character ``_``.
The name of each constant has to begin with a letter. The following names are used
by ImpactX, and cannot be used as user-defined constants: ``x``, ``y``, ``z``, ``t``, ``X``, ``Y``, ``Z``, ``T``.
The values of the constants can include the predefined ImpactX constants listed above as well as other user-defined constants.
For example:

* ``my_constants.my_alpha = 3.0``
* ``my_constants.my_beta = 12.e-6``
* ``my_constants.abc = 1.23e10``


Coordinates
^^^^^^^^^^^

Besides, for profiles that depend on spatial coordinates, the parser will interpret some variables as spatial coordinates.
These are specified in the input parameter, i.e., ``field_function(x,y,z)`` or ``field_function(X,Y,T)``.

The parser reads python-style expressions between double quotes, for instance
``"a0*x**2 * (1-y*1.e2) * (x>0)"`` is a valid expression where ``a0`` is a
user-defined constant (see above) and ``x`` and ``y`` are spatial coordinates. The names are case sensitive. The factor
``(x>0)`` is ``1`` where ``x>0`` and ``0`` where ``x<=0``. It allows the user to
define functions by intervals.
Alternatively the expression above can be written as ``if(x>0, a0*x**2 * (1-y*1.e2), 0)``.


.. _running-cpp-parameters-diagnostics:

Diagnostics and output
----------------------

.. pp:param:: diag.enable
    :type: ``boolean``
    :optional:
    :default: ``true``

    Enable or disable diagnostics generally.
    Disabling this is mostly used for benchmarking.

    This option is ignored for the openPMD output elements (remove them from the lattice to disable).

.. pp:param:: diag.slice_step_diagnostics
    :type: ``boolean``
    :optional:
    :default: ``false``

    By default, diagnostics are computed and written at the beginning and end of the simulation.
    Enabling this flag will write diagnostics at every step and slice step.

.. pp:param:: diag.file_prefix
    :type: ``string``
    :optional:
    :default: ``diags``

    Root directory for diagnostic output.
    By default, diagnostics are written in the folder ``diags/``.

    Set to an empty string or ``.`` to write diagnostics in the current working directory.

    If a directory at :pp:param:`diag.file_prefix` already exists when a simulation starts,
    ImpactX renames it to ``<diag.file_prefix>.old.<suffix>`` to preserve prior results.
    This is skipped when :pp:param:`diag.file_prefix` resolves to the current working directory,
    the root directory, or an ancestor of the current working directory; in those cases
    new output is written alongside existing files.

.. pp:param:: diag.file_min_digits
    :type: ``integer``
    :optional:
    :default: ``6``

    The minimum number of digits used for the step number appended to the diagnostic file names.

.. pp:param:: diag.backend
    :type: ``string``
    :default: ``default``

    Diagnostics for particles lost in apertures, stored as ``<diag.file_prefix>/openPMD/particles_lost.*`` at the end of the simulation.
    See the ``beam_monitor`` element for backend values.

.. pp:param:: diag.eigenemittances
    :type: ``boolean``
    :optional:
    :default: ``false``

    If this flag is enabled, the 3 eigenemittances of the 6D beam distribution are computed and written as diagnostics.
    This flag is disabled by default to reduce computational cost.


.. _running-cpp-parameters-cp-restart:

Checkpoints and restart
-----------------------

.. note::

   Future version of ImpactX will support checkpoints/restart.
   This is not yet implemented.


Intervals parser
----------------

.. note::

   TODO :-)

ImpactX can parse time step interval expressions of the form ``start:stop:period``, e.g.
``1:2:3, 4::, 5:6, :, ::10``.
A comma is used as a separator between groups of intervals, which we call slices.
The resulting time steps are the `union set <https://en.wikipedia.org/wiki/Union_(set_theory)>`_ of all given slices.
White spaces are ignored.
A single slice can have 0, 1 or 2 colons ``:``, just as `numpy slices <https://numpy.org/doc/stable/reference/generated/numpy.s_.html>`_, but with inclusive upper bound for ``stop``.

* For 0 colon the given value is the period

* For 1 colon the given string is of the type ``start:stop``

* For 2 colons the given string is of the type ``start:stop:period``

Any value that is not given is set to default.
Default is ``0`` for the start, ``std::numeric_limits<int>::max()`` for the stop and ``1`` for the
period.
For the 1 and 2 colon syntax, actually having values in the string is optional
(this means that ``::5``, ``100 ::10`` and ``100 :`` are all valid syntaxes).

All values can be expressions that will be parsed in the same way as other integer input parameters.

**Examples**

* ``something_intervals = 50`` -> do something at timesteps 0, 50, 100, 150, etc.
  (equivalent to ``something_intervals = ::50``)

* ``something_intervals = 300:600:100`` -> do something at timesteps 300, 400, 500 and 600.

* ``something_intervals = 300::50`` -> do something at timesteps 300, 350, 400, 450, etc.

* ``something_intervals = 105:108,205:208`` -> do something at timesteps 105, 106, 107, 108,
  205, 206, 207 and 208. (equivalent to ``something_intervals = 105 : 108 : , 205 : 208 :``)

* ``something_intervals = :`` or  ``something_intervals = ::`` -> do something at every timestep.

* ``something_intervals = 167:167,253:253,275:425:50`` do something at timesteps 167, 253, 275,
  325, 375 and 425.

This is essentially the python slicing syntax except that the stop is inclusive
(``0:100`` contains 100) and that no colon means that the given value is the period.

Note that if a given period is zero or negative, the corresponding slice is disregarded.
For example, ``something_intervals = -1`` deactivates ``something`` and
``something_intervals = ::-1,100:1000:25`` is equivalent to ``something_intervals = 100:1000:25``.


.. _running-cpp-parameters-overall:

Overall simulation parameters
-----------------------------

.. pp:param:: amrex.abort_on_out_of_gpu_memory
    :type: ``0`` or ``1``
    :default: ``1`` (true)

    When running on GPUs, memory that does not fit on the device will be automatically swapped to host memory when this option is set to ``0``.
    This will cause severe performance drops.
    Note that even with this set to ``1`` ImpactX will not catch all out-of-memory events yet when operating close to maximum device memory.
    `Please also see the documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters>`__.

.. pp:param:: amrex.the_arena_is_managed
    :type: ``0`` or ``1``
    :default: ``0`` (false)

    When running on GPUs, device memory that is accessed from the host will automatically be transferred with managed memory.
    This is useful for convenience during development, but has sometimes severe performance and memory footprint implications if relied on (and sometimes vendor bugs).
    For all regular ImpactX operations, we therefore do explicit memory transfers without the need for managed memory.
    `Please also see the documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters>`__.

.. pp:param:: amrex.omp_threads
    :type: ``system``, ``nosmt`` or positive integer
    :default: ``nosmt``

    An integer number can be set in lieu of the ``OMP_NUM_THREADS`` environment variable to control the number of OpenMP threads to use for the ``OMP`` compute backend on CPUs.
    By default, we use the ``nosmt`` option, which overwrites the OpenMP default of spawning one thread per logical CPU core, and instead only spawns a number of threads equal to the number of physical CPU cores on the machine.
    If set, the environment variable ``OMP_NUM_THREADS`` takes precedence over ``system`` and ``nosmt``, but not over integer numbers set in this option.

.. pp:param:: amrex.abort_on_unused_inputs
    :type: ``0`` or ``1``
    :default: ``0`` (false)

    When set to ``1``, this option causes the simulation to fail *after* its completion if there were unused parameters.
    It is mainly intended for continuous integration and automated testing to check that all tests and inputs are adapted to API changes.

.. pp:param:: impactx.always_warn_immediately
    :type: ``0`` or ``1``
    :default: ``0`` (false)

    If set to ``1``, ImpactX immediately prints every warning message as soon as it is generated.
    It is mainly intended for debug purposes, in case a simulation crashes before a global warning report can be printed.

.. pp:param:: impactx.abort_on_warning_threshold
    :type: ``string``
    :optional:

    Allowed values: ``low``, ``medium`` or ``high``.
    Optional threshold to abort as soon as a warning is raised.
    If the threshold is set, warning messages with priority greater than or equal to the threshold trigger an immediate abort.
    It is mainly intended for debug purposes, and is best used with :pp:param:`impactx.always_warn_immediately = 1`.
    For more information on the warning logger, see `this section <https://warpx.readthedocs.io/en/latest/developers/warning_logger.html>`__ of the WarpX documentation.

.. pp:param:: impactx.verbose
    :type: ``integer``
    :optional:
    :default: ``1``

    Controls how much information is printed to the terminal, when running ImpactX.
    Use ``0`` for silent; higher is more verbose.
