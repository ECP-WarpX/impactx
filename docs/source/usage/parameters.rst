.. _running-cpp-parameters:

Parameters: Inputs File
=======================

This documents on how to use ImpactX with an inputs file (``impactx input_file.in``).

.. note::
   The AMReX parser (see :ref:`running-cpp-parameters-parser`) is used for the right-hand-side of all input parameters that consist of one or more integers or floats, so expressions like ``<species_name>.density_max = "2.+1."`` and/or using user-defined constants are accepted.

.. _running-cpp-parameters-overall:

Overall simulation parameters
-----------------------------

* ``max_step`` (``integer``)
    The number of PIC cycles to perform.

* ``stop_time`` (``float``; in seconds)
    The maximum physical time of the simulation. Can be provided instead of ``max_step``. If both
    ``max_step`` and ``stop_time`` are provided, both criteria are used and the simulation stops
    when the first criterion is hit.

* ``amrex.abort_on_out_of_gpu_memory``  (``0`` or ``1``; default is ``1`` for true)
    When running on GPUs, memory that does not fit on the device will be automatically swapped to host memory when this option is set to ``0``.
    This will cause severe performance drops.
    Note that even with this set to ``1`` ImpactX will not catch all out-of-memory events yet when operating close to maximum device memory.
    `Please also see the documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters>`__.

* ``amrex.the_arena_is_managed``  (``0`` or ``1``; default is ``0`` for false)
    When running on GPUs, device memory that is accessed from the host will automatically be transferred with managed memory.
    This is useful for convenience during development, but has sometimes severe performance and memory footprint implications if relied on (and sometimes vendor bugs).
    For all regular ImpactX operations, we therefore do explicit memory transfers without the need for managed memory and thus changed the AMReX default to false.
    `Please also see the documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters>`__.

* ``amrex.omp_threads``  (``system``, ``nosmt`` or positive integer; default is ``nosmt``)
    An integer number can be set in lieu of the ``OMP_NUM_THREADS`` environment variable to control the number of OpenMP threads to use for the ``OMP`` compute backend on CPUs.
    By default, we use the ``nosmt`` option, which overwrites the OpenMP default of spawning one thread per logical CPU core, and instead only spawns a number of threads equal to the number of physical CPU cores on the machine.
    If set, the environment variable ``OMP_NUM_THREADS`` takes precedence over ``system`` and ``nosmt``, but not over integer numbers set in this option.

* ``amrex.abort_on_unused_inputs`` (``0`` or ``1``; default is ``0`` for false)
    When set to ``1``, this option causes the simulation to fail *after* its completion if there were unused parameters.
    It is mainly intended for continuous integration and automated testing to check that all tests and inputs are adapted to API changes.

* ``impactx.always_warn_immediately`` (``0`` or ``1``; default is ``0`` for false)
    If set to ``1``, ImpactX immediately prints every warning message as soon as it is generated.
    It is mainly intended for debug purposes, in case a simulation crashes before a global warning report can be printed.

* ``impactx.abort_on_warning_threshold`` (string: ``low``, ``medium`` or ``high``) optional
    Optional threshold to abort as soon as a warning is raised.
    If the threshold is set, warning messages with priority greater than or equal to the threshold trigger an immediate abort.
    It is mainly intended for debug purposes, and is best used with ``impactx.always_warn_immediately=1``.
    For more information on the warning logger, see `this section <https://warpx.readthedocs.io/en/latest/developers/warning_logger.html>`__ of the WarpX documentation.

.. _running-cpp-parameters-box:


Setting up the field mesh
-------------------------

ImpactX uses an AMReX grid of boxes to organize and parallelize the simulation domain.
These boxes also contain a field mesh, if space charge calculations are enabled.

* ``amr.n_cell`` (3 integers) optional (default: 1 `blocking_factor <https://amrex-codes.github.io/amrex/docs_html/GridCreation.html>`__ per MPI process)
    The number of grid points along each direction (on the **coarsest level**)

* ``amr.max_level`` (``integer``, default: ``0``)
    When using mesh refinement, the number of refinement levels that will be used.

    Use ``0`` in order to disable mesh refinement.

* ``amr.ref_ratio`` (``integer`` per refined level, default: ``2``)
    When using mesh refinement, this is the refinement ratio per level.
    With this option, all directions are fined by the same ratio.

* ``amr.ref_ratio_vect`` (3 integers for x,y,z per refined level)
    When using mesh refinement, this can be used to set the refinement ratio per direction and level, relative to the previous level.

    Example: for three levels, a value of ``2 2 4 8 8 16`` refines the first level by 2-fold in x and y and 4-fold in z compared to the coarsest level (level 0/mother grid); compared to the first level, the second level is refined 8-fold in x and y and 16-fold in z.

.. note::

   Field boundaries for space charge calculation are located at the outer ends of the field mesh.
   We currently assume `Dirichlet boundary conditions <https://en.wikipedia.org/wiki/Dirichlet_boundary_condition>`__ with zero potential (a mirror charge).
   Thus, to emulate open boundaries, consider adding enough vacuum padding to the beam.
   This will be improved in future versions.

.. note::

   Particles that move outside the simulation domain are removed.

* ``geometry.dynamic_size`` (``boolean``) optional (default: ``true`` for dynamic)
    Use dynamic (``true``) resizing of the field mesh, via ``geometry.prob_relative``, or static sizing (``false``), via ``geometry.prob_lo``/``geometry.prob_hi``.

* ``geometry.prob_relative`` (positive ``float`` array with ``amr.max_level`` entries, unitless) optional (default: ``3.0 1.0 1.0 ...``)
    By default, we dynamically extract the minimum and maximum of the particle positions in the beam.
    The field mesh spans, per direction, multiple times the maximum physical extent of beam particles, as given by this factor.
    The beam minimum and maximum extent are symmetrically padded by the mesh.
    For instance, ``1.2`` means the mesh will span 10% above and 10% below the beam;
    ``1.0`` means the beam is exactly covered with the mesh.

* ``geometry.prob_lo`` and ``geometry.prob_hi`` (3 floats, in meters) optional (required if ``geometry.dynamic_size`` is ``false``)
    The extent of the full simulation domain relative to the reference particle position.
    This can be used to explicitly size the simulation box and ignore ``geometry.prob_relative``.

    This box is rectangular, and thus its extent is given here by the coordinates of the lower corner (``geometry.prob_lo``) and upper corner (``geometry.prob_hi``).
    The first axis of the coordinates is x and the last is z.


.. _running-cpp-parameters-bc:

Domain Boundary Conditions
--------------------------

.. note::

   TODO :-)


.. _running-cpp-parameters-particle:

Initial Beam Distributions
--------------------------

* ``beam.npart`` (``integer``)
  number of weighted simulation particles

* ``beam.units`` (``string``)
  currently, only ``static`` is supported.

* ``beam.kin_energy`` (``float``, in MeV)
  beam kinetic energy

* ``beam.charge`` (``float``, in C)
  bunch charge

* ``beam.particle`` (``string``)
  particle type: currently either ``electron``, ``positron`` or ``proton``

* ``beam.distribution`` (``string``)
    Indicates the initial distribution type.
    This should be one of:

    * ``waterbag`` for initial Waterbag distribution.
      With additional parameters:

        * ``beam.sigmaX`` (``float``, in meters) rms X
        * ``beam.sigmaY`` (``float``, in meters) rms Y
        * ``beam.sigmaT`` (``float``, in radian) rms normalized time difference T
        * ``beam.sigmaPx`` (``float``, in momentum) rms Px
        * ``beam.sigmaPy`` (``float``, in momentum) rms Py
        * ``beam.sigmaPt`` (``float``, in energy deviation) rms Pt
        * ``beam.muxpx`` (``float``, dimensionless, default: ``0``) correlation X-Px
        * ``beam.muypy`` (``float``, dimensionless, default: ``0``) correlation Y-Py
        * ``beam.mutpt`` (``float``, dimensionless, default: ``0``) correlation T-Pt

    * ``kurth6d`` for initial 6D Kurth distribution.
      With additional parameters:

        * ``beam.sigmaX`` (``float``, in meters) rms X
        * ``beam.sigmaY`` (``float``, in meters) rms Y
        * ``beam.sigmaT`` (``float``, in radian) rms normalized time difference T
        * ``beam.sigmaPx`` (``float``, in momentum) rms Px
        * ``beam.sigmaPy`` (``float``, in momentum) rms Py
        * ``beam.sigmaPt`` (``float``, in energy deviation) rms Pt
        * ``beam.muxpx`` (``float``, dimensionless, default: ``0``) correlation X-Px
        * ``beam.muypy`` (``float``, dimensionless, default: ``0``) correlation Y-Py
        * ``beam.mutpt`` (``float``, dimensionless, default: ``0``) correlation T-Pt

    * ``gaussian`` for initial 6D Gaussian (normal) distribution.
      With additional parameters:

        * ``beam.sigmaX`` (``float``, in meters) rms X
        * ``beam.sigmaY`` (``float``, in meters) rms Y
        * ``beam.sigmaT`` (``float``, in radian) rms normalized time difference T
        * ``beam.sigmaPx`` (``float``, in momentum) rms Px
        * ``beam.sigmaPy`` (``float``, in momentum) rms Py
        * ``beam.sigmaPt`` (``float``, in energy deviation) rms Pt
        * ``beam.muxpx`` (``float``, dimensionless, default: ``0``) correlation X-Px
        * ``beam.muypy`` (``float``, dimensionless, default: ``0``) correlation Y-Py
        * ``beam.mutpt`` (``float``, dimensionless, default: ``0``) correlation T-Pt

    * ``kvdist`` for initial K-V distribution in the transverse plane.
      The distribution is uniform in t and Gaussian in pt.
      With additional parameters:

        * ``beam.sigmaX`` (``float``, in meters) rms X
        * ``beam.sigmaY`` (``float``, in meters) rms Y
        * ``beam.sigmaT`` (``float``, in radian) rms normalized time difference T
        * ``beam.sigmaPx`` (``float``, in momentum) rms Px
        * ``beam.sigmaPy`` (``float``, in momentum) rms Py
        * ``beam.sigmaPt`` (``float``, in energy deviation) rms Pt
        * ``beam.muxpx`` (``float``, dimensionless, default: ``0``) correlation X-Px
        * ``beam.muypy`` (``float``, dimensionless, default: ``0``) correlation Y-Py
        * ``beam.mutpt`` (``float``, dimensionless, default: ``0``) correlation T-Pt

    * ``kurth4d`` for initial 4D Kurth distribution in the transverse plane.
      The distribution is uniform in t and Gaussian in pt.
      With additional parameters:

        * ``beam.sigmaX`` (``float``, in meters) rms X
        * ``beam.sigmaY`` (``float``, in meters) rms Y
        * ``beam.sigmaT`` (``float``, in radian) rms normalized time difference T
        * ``beam.sigmaPx`` (``float``, in momentum) rms Px
        * ``beam.sigmaPy`` (``float``, in momentum) rms Py
        * ``beam.sigmaPt`` (``float``, in energy deviation) rms Pt
        * ``beam.muxpx`` (``float``, dimensionless, default: ``0``) correlation X-Px
        * ``beam.muypy`` (``float``, dimensionless, default: ``0``) correlation Y-Py
        * ``beam.mutpt`` (``float``, dimensionless, default: ``0``) correlation T-Pt

    * ``semigaussian`` for initial Semi-Gaussian distribution.  The distribution is uniform within a cylinder in (x,y,z) and Gaussian in momenta (px,py,pt).
      With additional parameters:

        * ``beam.sigmaX`` (``float``, in meters) rms X
        * ``beam.sigmaY`` (``float``, in meters) rms Y
        * ``beam.sigmaT`` (``float``, in radian) rms normalized time difference T
        * ``beam.sigmaPx`` (``float``, in momentum) rms Px
        * ``beam.sigmaPy`` (``float``, in momentum) rms Py
        * ``beam.sigmaPt`` (``float``, in energy deviation) rms Pt
        * ``beam.muxpx`` (``float``, dimensionless, default: ``0``) correlation X-Px
        * ``beam.muypy`` (``float``, dimensionless, default: ``0``) correlation Y-Py
        * ``beam.mutpt`` (``float``, dimensionless, default: ``0``) correlation T-Pt

    * ``triangle`` a triangle distribution for laser-plasma acceleration related applications.
      A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).
      The transverse distribution is a 4D waterbag.
      With additional parameters:

        * ``beam.sigmaX`` (``float``, in meters) rms X
        * ``beam.sigmaY`` (``float``, in meters) rms Y
        * ``beam.sigmaT`` (``float``, in radian) rms normalized time difference T
        * ``beam.sigmaPx`` (``float``, in momentum) rms Px
        * ``beam.sigmaPy`` (``float``, in momentum) rms Py
        * ``beam.sigmaPt`` (``float``, in energy deviation) rms Pt
        * ``beam.muxpx`` (``float``, dimensionless, default: ``0``) correlation X-Px
        * ``beam.muypy`` (``float``, dimensionless, default: ``0``) correlation Y-Py
        * ``beam.mutpt`` (``float``, dimensionless, default: ``0``) correlation T-Pt

    * ``thermal`` for a 6D stationary thermal or bithermal distribution.
      This distribution type is described, for example in:
      R. D. Ryne et al, "A Test Suite of Space-Charge Problems for Code Benchmarking", in Proc. EPAC2004, Lucerne, Switzerland.
      C. E. Mitchell et al, "ImpactX Modeling of Benchmark Tests for Space Charge Validation", in Proc. HB2023, Geneva, Switzerland.
      With additional parameters:

        * ``beam.k`` (``float``, in inverse meters) external focusing strength
        * ``beam.kT`` (``float``, dimensionless) temperature of core population
           = < p_x^2 > = < p_y^2 >, where all momenta are normalized by the reference momentum
        * ``beam.kT_halo`` (``float``, dimensionless, default ``kT``) temperature of halo population
        * ``beam.normalize`` (``float``, dimensionless) normalizing constant for core population
        * ``beam.normalize_halo`` (``float``, dimensionless) normalizing constant for halo population
        * ``beam.halo`` (``float``, dimensionless) fraction of charge in halo


.. _running-cpp-parameters-lattice:

Lattice Elements
----------------

* ``lattice.elements`` (``list of strings``) optional (default: no elements)
    A list of names (one name per lattice element), in the order that they appear in the lattice.

* ``lattice.periods`` (``integer``) optional (default: ``1``)
    The number of periods to repeat the lattice.

* ``lattice.reverse`` (``boolean``) optional (default: ``false``)
    Reverse the list of elements in the lattice.
    If ``reverse`` and ``periods`` both appear, then ``reverse`` is applied before ``periods``.

* ``lattice.nslice`` (``integer``) optional (default: ``1``)
    A positive integer specifying the number of slices used for the application of
    space charge in all elements; overwritten by element parameter "nslice"

* ``<element_name>.type`` (``string``)
    Indicates the element type for this lattice element. This should be one of:

         * ``cfbend`` for a combined function bending magnet. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.rc`` (``float``, in meters) the bend radius
            * ``<element_name>.k`` (``float``, in inverse meters squared) the quadrupole strength
                = (magnetic field gradient in T/m) / (magnetic rigidity in T-m)

              * k > 0 horizontal focusing
              * k < 0 horizontal defocusing

            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``drift`` for a free drift. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``drift_chromatic`` for a free drift, with chromatic effects included.
           The Hamiltonian is expanded through second order in the transverse variables (x,px,y,py), with the exact pt dependence retained.
           This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``drift_exact`` for a free drift, using the exact nonlinear map. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``quad`` for a quadrupole. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.k`` (``float``, in inverse meters squared) the quadrupole strength
                = (magnetic field gradient in T/m) / (magnetic rigidity in T-m)

              * k > 0 horizontal focusing
              * k < 0 horizontal defocusing

            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``quad_chromatic`` for A Quadrupole magnet, with chromatic effects included.
           The Hamiltonian is expanded through second order in the transverse variables (x,px,y,py), with the exact pt dependence retained.
           This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.k`` (``float``, in inverse meters squared OR in T/m) the quadrupole strength
                = (magnetic field gradient in T/m) / (magnetic rigidity in T-m) - if units = 0

             OR = magnetic field gradient in T/m - if units = 1

              * k > 0 horizontal focusing
              * k < 0 horizontal defocusing

            * ``<element_name>.units`` (``integer``) specification of units (default: ``0``)
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``quadrupole_softedge`` for a soft-edge quadrupole. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.gscale`` (``float``, in inverse meters) Scaling factor for on-axis magnetic field gradient
            * ``<element_name>.cos_coefficients`` (array of ``float``) cos coefficients in Fourier expansion of the on-axis field gradient
              (optional); default is a tanh fringe field model from `MaryLie 3.0 <http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__
            * ``<element_name>.sin_coefficients`` (array of ``float``) sin coefficients in Fourier expansion of the on-axis field gradient
              (optional); default is a tanh fringe field model from `MaryLie 3.0 <http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.mapsteps`` (``integer``) number of integration steps per slice used for map and reference particle push in applied fields
               (default: ``1``)
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``sbend`` for a bending magnet. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.rc`` (``float``, in meters) the bend radius
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``sbend_exact`` for a bending magnet using the exact nonlinear map for the bend body. The map corresponds to the map described in:
            D. L. Bruhwiler et al, in Proc. of EPAC 98, pp. 1171-1173 (1998), E. Forest et al, Part. Accel. 45, pp. 65-94 (1994).  The model
            consists of a uniform bending field B_y with a hard edge.  Pole faces are normal to the entry and exit velocity of the reference
            particle.  This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.phi`` (``float``, in degrees) the bend angle
            * ``<element_name>.B`` (``float``, in Tesla) the bend magnetic field; when B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``solenoid`` for an ideal hard-edge solenoid magnet. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.ks`` (``float``, in meters) Solenoid strength in m^(-1) (MADX convention)
                  = (magnetic field Bz in T) / (rigidity in T-m)
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``solenoid_softedge`` for a soft-edge solenoid. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.bscale`` (``float``, in inverse meters) Scaling factor for on-axis magnetic field Bz
            * ``<element_name>.cos_coefficients`` (array of ``float``) cos coefficients in Fourier expansion of the on-axis magnetic field Bz
              (optional); default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__
            * ``<element_name>.sin_coefficients`` (array of ``float``) sin coefficients in Fourier expansion of the on-axis magnetic field Bz
              (optional); default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.mapsteps`` (``integer``) number of integration steps per slice used for map and reference particle push in applied fields (default: ``1``)
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``dipedge`` for dipole edge focusing. This requires these additional parameters:

            * ``<element_name>.psi`` (``float``, in radians) the pole face rotation angle
            * ``<element_name>.rc`` (``float``, in meters) the bend radius
            * ``<element_name>.g`` (``float``, in meters) the gap size
            * ``<element_name>.K2`` (``float``, dimensionless) normalized field integral for fringe field
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``constf`` for a constant focusing element. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.kx`` (``float``, in 1/meters) the horizontal focusing strength
            * ``<element_name>.ky`` (``float``, in 1/meters) the vertical focusing strength
            * ``<element_name>.kt`` (``float``, in 1/meters) the longitudinal focusing strength
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``rfcavity`` a radiofrequency cavity.
          This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.escale`` (``float``, in 1/m) scaling factor for on-axis RF electric field
                = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV)
            * ``<element_name>.freq`` (``float``, in Hz) RF frequency
            * ``<element_name>.phase`` (``float``, in degrees) RF driven phase
            * ``<element_name>.cos_coefficients`` (array of ``float``) cosine coefficients in Fourier expansion of on-axis electric field Ez (optional); default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__
            * ``<element_name>.cos_coefficients`` (array of ``float``) sine coefficients in Fourier expansion of on-axis electric field Ez (optional); default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.mapsteps`` (``integer``) number of integration steps per slice used for map and reference particle push in applied fields (default: ``1``)
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``buncher`` for a short RF cavity (linear) bunching element.
          This requires these additional parameters:

            * ``<element_name>.V`` (``float``, dimensionless) normalized voltage drop across the cavity
                = (maximum voltage drop in Volts) / (speed of light in m/s * magnetic rigidity in T-m)
            * ``<element_name>.k`` (``float``, in 1/meters) the RF wavenumber
                = 2*pi/(RF wavelength in m)
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``shortrf`` for a short RF cavity element.
          This requires these additional parameters:

            * ``<element_name>.V`` (``float``, dimensionless) normalized voltage drop across the cavity
                = (maximum energy gain in MeV) / (particle rest energy in MeV)
            * ``<element_name>.freq`` (``float``, in Hz) the RF frequency
            * ``<element_name>.phase`` (``float``, in degrees) the synchronous RF phase

                phase = 0: maximum energy gain (on-crest)

                phase = -90 deg:  zero energy gain for bunching

                phase = 90 deg:  zero energy gain for debunching
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``uniform_acc_chromatic`` for a region of uniform acceleration, with chromatic effects included.
           The Hamiltonian is expanded through second order in the transverse variables (x,px,y,py), with the exact pt dependence retained.
           This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length
            * ``<element_name>.ez`` (``float``, in inverse meters) the electric field strength
                = (particle charge in C * electric field Ez in V/m) / (particle mass in kg * (speed of light in m/s)^2)
            * ``<element_name>.bz`` (``float``, in inverse meters) the magnetic field strength
                = (particle charge in C * magnetic field Bz in T) / (particle mass in kg * speed of light in m/s)
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane
            * ``<element_name>.nslice`` (``integer``) number of slices used for the application of space charge (default: ``1``)

        * ``multipole`` for a thin multipole element.
          This requires these additional parameters:

            * ``<element_name>.multipole`` (``integer``, dimensionless) order of multipole
                (m = 1) dipole, (m = 2) quadrupole, (m = 3) sextupole, etc.

            * ``<element_name>.k_normal`` (``float``, in 1/meters^m) integrated normal multipole coefficient (MAD-X convention)
                = 1/(magnetic rigidity in T-m) * (derivative of order m-1 of By with respect to x)
            * ``<element_name>.k_skew`` (``float``, in 1/meters^m) integrated skew multipole strength (MAD-X convention)
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``nonlinear_lens`` for a thin IOTA nonlinear lens element.
          This requires these additional parameters:

            * ``<element_name>.knll`` (``float``, in meters) integrated strength of the lens segment (MAD-X convention)
                = dimensionless lens strength * c parameter**2 * length / Twiss beta
            * ``<element_name>.cnll`` (``float``, in meters) distance of the singularities from the origin (MAD-X convention)
                = c parameter * sqrt(Twiss beta)
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``prot`` for an exact pole-face rotation in the x-z plane. This requires these additional parameters:

            * ``<element_name>.phi_in`` (``float``, in degrees) angle of the reference particle with respect to the longitudinal (z) axis in the original frame
            * ``<element_name>.phi_out`` (``float``, in degrees) angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame

        * ``kicker`` for a thin transverse kicker. This requires these additional parameters:

            * ``<element_name>.xkick`` (``float``, dimensionless OR in T-m) the horizontal kick strength
            * ``<element_name>.ykick`` (``float``, dimensionless OR in T-m) the vertical kick strength
            * ``<element_name>.units`` (``string``) specification of units: ``dimensionless`` (default, in units of the magnetic rigidity of the reference particle) or ``T-m``
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``thin_dipole`` for a thin dipole element.
          This requires these additional parameters:

            * ``<element_name>.theta`` (``float``, in degrees) dipole bend angle
            * ``<element_name>.rc`` (``float``, in meters) effective radius of curvature
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``aperture`` for a thin collimator element applying a transverse aperture boundary.
          This requires these additional parameters:

            * ``<element_name>.xmax`` (``float``, in meters) maximum value of the horizontal coordinate
            * ``<element_name>.ymax`` (``float``, in meters) maximum value of the vertical coordinate
            * ``<element_name>.shape`` (``string``) shape of the aperture boundary: ``rectangular`` (default) or ``elliptical``
            * ``<element_name>.dx`` (``float``, in meters) horizontal translation error
            * ``<element_name>.dy`` (``float``, in meters) vertical translation error
            * ``<element_name>.rotation`` (``float``, in degrees) rotation error in the transverse plane

        * ``beam_monitor`` a beam monitor, writing all beam particles at fixed ``s`` to openPMD files.
          If the same element name is used multiple times, then an output series is created with multiple outputs.

            * ``<element_name>.name`` (``string``, default value: ``<element_name>``)

                The output series name to use.
                By default, output is created under ``diags/openPMD/<element_name>.<backend>``.

            * ``<element_name>.backend`` (``string``, default value: ``default``)

                `I/O backend <https://openpmd-api.readthedocs.io/en/latest/backends/overview.html>`_ for `openPMD <https://www.openPMD.org>`_ data dumps.
                ``bp`` is the `ADIOS2 I/O library <https://csmd.ornl.gov/adios>`_, ``h5`` is the `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_, and ``json`` is a `simple text format <https://en.wikipedia.org/wiki/JSON>`_.
                ``json`` only works with serial/single-rank jobs.
                By default, the first available backend in the order given above is taken.

            * ``<element_name>.encoding`` (``string``, default value: ``g``)

                openPMD `iteration encoding <https://openpmd-api.readthedocs.io/en/0.14.0/usage/concepts.html#iteration-and-series>`__: (v)ariable based, (f)ile based, (g)roup based (default)
                variable based is an `experimental feature with ADIOS2 <https://openpmd-api.readthedocs.io/en/0.14.0/backends/adios2.html#experimental-new-adios2-schema>`__.

        * ``line`` a sub-lattice (line) of elements to append to the lattice.

            * ``<element_name>.elements`` (``list of strings``) optional (default: no elements)
              A list of names (one name per lattice element), in the order that they appear in the lattice.

            * ``<element_name>.reverse`` (``boolean``) optional (default: ``false``)
              Reverse the list of elements in the line before appending to the lattice.

            * ``<element_name>.repeat`` (``integer``) optional (default: ``1``)
              Repeat the line multiple times before appending to the lattice.
              Note: If ``reverse`` and ``repeat`` both appear, then ``reverse`` is applied before ``repeat``.


.. _running-cpp-parameters-parallelization:

Distribution across MPI ranks and parallelization
-------------------------------------------------

* ``amr.max_grid_size`` (``integer``) optional (default: ``128``)
    Maximum allowable size of each **subdomain**
    (expressed in number of grid points, in each direction).
    Each subdomain has its own ghost cells, and can be handled by a
    different MPI rank ; several OpenMP threads can work simultaneously on the
    same subdomain.

    If ``max_grid_size`` is such that the total number of subdomains is
    **larger** that the number of MPI ranks used, than some MPI ranks
    will handle several subdomains, thereby providing additional flexibility
    for **load balancing**.

    When using mesh refinement, this number applies to the subdomains
    of the coarsest level, but also to any of the finer level.


.. _running-cpp-parameters-parser:

Math parser and user-defined constants
--------------------------------------

ImpactX uses AMReX's math parser that reads expressions in the input file.
It can be used in all input parameters that consist of one or more integers or floats.
Integer input expecting boolean, 0 or 1, are not parsed.
Note that when multiple values are expected, the expressions are space delimited.
For integer input values, the expressions are evaluated as real numbers and the final result rounded to the nearest integer.
See `this section <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parser>`_ of the AMReX documentation for a complete list of functions supported by the math parser.

ImpactX constants
^^^^^^^^^^^^^^^^^

ImpactX will provide a few pre-defined constants, that can be used for any parameter that consists of one or more floats.

.. note::

   Develop, such as:

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

   See in WarpX the file ``Source/Utils/WarpXConst.H`` for the values.

User-defined constants
^^^^^^^^^^^^^^^^^^^^^^

Users can define their own constants in the input file.
These constants can be used for any parameter that consists of one or more integers or floats.
User-defined constant names can contain only letters, numbers and the character ``_``.
The name of each constant has to begin with a letter. The following names are used
by ImpactX, and cannot be used as user-defined constants: ``x``, ``y``, ``z``, ``X``, ``Y``, ``t``.
The values of the constants can include the predefined ImpactX constants listed above as well as other user-defined constants.
For example:

* ``my_constants.a0 = 3.0``
* ``my_constants.z_plateau = 150.e-6``
* ``my_constants.n0 = 1.e22``
* ``my_constants.wp = sqrt(n0*q_e**2/(epsilon0*m_e))``

Coordinates
^^^^^^^^^^^

Besides, for profiles that depend on spatial coordinates (the plasma momentum distribution or the laser field, see below ``Particle initialization`` and ``Laser initialization``), the parser will interpret some variables as spatial coordinates.
These are specified in the input parameter, i.e., ``density_function(x,y,z)`` and ``field_function(X,Y,t)``.

The parser reads python-style expressions between double quotes, for instance
``"a0*x**2 * (1-y*1.e2) * (x>0)"`` is a valid expression where ``a0`` is a
user-defined constant (see above) and ``x`` and ``y`` are spatial coordinates. The names are case sensitive. The factor
``(x>0)`` is ``1`` where ``x>0`` and ``0`` where ``x<=0``. It allows the user to
define functions by intervals.
Alternatively the expression above can be written as ``if(x>0, a0*x**2 * (1-y*1.e2), 0)``.


.. _running-cpp-parameters-numerics:

Numerics and algorithms
-----------------------

* ``algo.particle_shape`` (``integer``; ``1``, ``2``, or ``3``)
    The order of the shape factors (splines) for the macro-particles along all spatial directions: `1` for linear, `2` for quadratic, `3` for cubic.
    Low-order shape factors result in faster simulations, but may lead to more noisy results.
    High-order shape factors are computationally more expensive, but may increase the overall accuracy of the results.
    For production runs it is generally safer to use high-order shape factors, such as cubic order.

* ``algo.space_charge`` (``boolean``, optional, default: ``false``)
    Whether to calculate space charge effects.

* ``algo.mlmg_relative_tolerance`` (``float``, optional, default: ``1.e-7``)
    The relative precision with which the electrostatic space-charge fields should be calculated.
    More specifically, the space-charge fields are computed with an iterative Multi-Level Multi-Grid (MLMG) solver.
    This solver can fail to reach the default precision within a reasonable time.

* ``algo.mlmg_absolute_tolerance`` (``float``, optional, default: ``0``, which means: ignored)
    The absolute tolerance with which the space-charge fields should be calculated in units of V/m^2.
    More specifically, the acceptable residual with which the solution can be considered converged.
    In general this should be left as the default, but in cases where the simulation state changes very
    little between steps it can occur that the initial guess for the MLMG solver is so close to the
    converged value that it fails to improve that solution sufficiently to reach the
    mlmg_relative_tolerance value."

* ``algo.mlmg_max_iters`` (``integer``, optional, default: ``100``)
    Maximum number of iterations used for MLMG solver for space-charge fields calculation.
    In case if MLMG converges but fails to reach the desired self_fields_required_precision,
    this parameter may be increased.

* ``algo.mlmg_verbosity`` (``integer``, optional, default: ``1``)
    The verbosity used for MLMG solver for space-charge fields calculation.
    Currently MLMG solver looks for verbosity levels from 0-5.
    A higher number results in more verbose output.

.. _running-cpp-parameters-diagnostics:

Diagnostics and output
----------------------

* ``diag.enable`` (``boolean``, optional, default: ``true``)
  Enable or disable diagnostics generally.
  Disabling this is mostly used for benchmarking.

  This option is ignored for the openPMD output elements (remove them from the lattice to disable).

* ``diag.slice_step_diagnostics`` (``boolean``, optional, default: ``false``)
  By default, diagnostics is performed at the beginning and end of the simulation.
  Enabling this flag will write diagnostics every step and slice step

* ``diag.file_min_digits`` (``integer``, optional, default: ``6``)
    The minimum number of digits used for the step number appended to the diagnostic file names.

* ``diag.backend`` (``string``, default value: ``default``)

  Diagnostics for particles lost in apertures, stored as ``diags/openPMD/particles_lost.*`` at the end of the simulation.
  See the ``beam_monitor`` element for backend values.

.. _running-cpp-parameters-diagnostics-reduced:

Reduced Diagnostics
^^^^^^^^^^^^^^^^^^^

Reduced diagnostics allow the user to compute some reduced quantity (invariants of motion, particle temperature, max of a field, ...) and write a small amount of data to text files.
Reduced diagnostics are run *in situ* with the simulation.

Diagnostics related to integrable optics in the IOTA nonlinear magnetic insert element:

* ``diag.alpha`` (``float``, unitless) Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
  Horizontal and vertical values must be equal.

* ``diag.beta`` (``float``, meters) Twiss beta of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
  Horizontal and vertical values must be equal.

* ``diag.tn`` (``float``, unitless) dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I.

* ``diag.cn`` (``float``, meters^(1/2)) scale factor of the IOTA nonlinear magnetic insert element used for computing H and I.


.. _running-cpp-parameters-diagnostics-insitu:

In-situ visualization
^^^^^^^^^^^^^^^^^^^^^

.. note::

   TODO :-)

.. _running-cpp-parameters-diagnostics-full:

.. note::

   TODO :-)

.. _running-cpp-parameters-cp-restart:

Checkpoints and restart
-----------------------

.. note::

   ImpactX will support checkpoints/restart via AMReX.
   The checkpoint capability can be turned with regular diagnostics: ``<diag_name>.format = checkpoint``.

   * ``amr.restart`` (`string`)
       Name of the checkpoint file to restart from. Returns an error if the folder does not exist
       or if it is not properly formatted.

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
