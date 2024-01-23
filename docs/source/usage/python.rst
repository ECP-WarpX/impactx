.. _usage-picmi:

Parameters: Python
==================

This documents on how to use ImpactX as a Python script (``python3 run_script.py``).

General
-------

.. py:class:: impactx.ImpactX

   This is the central simulation class.

   .. py:property:: particle_shape

      Control the particle B-spline order.

      The order of the shape factors (splines) for the macro-particles along all spatial directions: `1` for linear, `2` for quadratic, `3` for cubic.
      Low-order shape factors result in faster simulations, but may lead to more noisy results.
      High-order shape factors are computationally more expensive, but may increase the overall accuracy of the results.
      For production runs it is generally safer to use high-order shape factors, such as cubic order.

      :param int order: B-spline order ``1``, ``2``, or ``3``

   .. py:property:: n_cell

      The number of grid points along each direction on the coarsest level.

   .. py:property:: max_level

      The maximum mesh-refinement level for the simulation.

   .. py:property:: finest_level

      The currently finest level of mesh-refinement used.
      This is always less or equal to :py:attr:`~max_level`.

   .. py:property:: domain

      The physical extent of the full simulation domain, relative to the reference particle position, in meters.
      When set, turns ``dynamic_size`` to ``False``.

      Note: particles that move outside the simulation domain are removed.

   .. py:property:: prob_relative

      This is a list with ``amr.max_level`` + 1 entries.

      By default, we dynamically extract the minimum and maximum of the particle positions in the beam.
      The field mesh spans, per direction, multiple times the maximum physical extent of beam particles, as given by this factor.
      The beam minimum and maximum extent are symmetrically padded by the mesh.
      For instance, ``1.2`` means the mesh will span 10% above and 10% below the beam;
      ``1.0`` means the beam is exactly covered with the mesh.

      Default: ``[3.0 1.0 1.0 ...]``.
      When set, turns ``dynamic_size`` to ``True``.

   .. py:property:: dynamic_size

      Use dynamic (``True``) resizing of the field mesh or static sizing (``False``).

   .. py:property:: space_charge

      Enable (``True``) or disable (``False``) space charge calculations (default: ``False``).

      Whether to calculate space charge effects.

   .. py:property:: mlmg_relative_tolerance

      Default: ``1.e-7``

      The relative precision with which the electrostatic space-charge fields should be calculated.
      More specifically, the space-charge fields are computed with an iterative Multi-Level Multi-Grid (MLMG) solver.
      This solver can fail to reach the default precision within a reasonable time.

   .. py:property:: mlmg_absolute_tolerance

      Default: ``0``, which means: ignored

      The absolute tolerance with which the space-charge fields should be calculated in units of :math:`V/m^2`.
      More specifically, the acceptable residual with which the solution can be considered converged.
      In general this should be left as the default, but in cases where the simulation state changes very
      little between steps it can occur that the initial guess for the MLMG solver is so close to the
      converged value that it fails to improve that solution sufficiently to reach the ``mlmg_relative_tolerance`` value.

   .. py:property:: mlmg_max_iters

      Default: ``100``

      Maximum number of iterations used for MLMG solver for space-charge fields calculation.
      In case if MLMG converges but fails to reach the desired self_fields_required_precision, this parameter may be increased.

   .. py:property:: mlmg_verbosity

      Default: ``1``

      The verbosity used for MLMG solver for space-charge fields calculation.
      Currently MLMG solver looks for verbosity levels from 0-5.
      A higher number results in more verbose output.

   .. py:property:: diagnostics

      Enable (``True``) or disable (``False``) diagnostics generally (default: ``True``).
      Disabling this is mostly used for benchmarking.

   .. py:property:: slice_step_diagnostics

      Enable (``True``) or disable (``False``) diagnostics every slice step in elements  (default: ``True``).

      By default, diagnostics is performed at the beginning and end of the simulation.
      Enabling this flag will write diagnostics every step and slice step.

   .. py:property:: diag_file_min_digits

      The minimum number of digits (default: ``6``) used for the step
      number appended to the diagnostic file names.

   .. py:method:: set_diag_iota_invariants(alpha, beta, tn, cn)

      Set the parameters of the IOTA nonlinear lens invariants diagnostics.

      :param float alpha: Twiss alpha
      :param float beta: Twiss beta (m)
      :param float tn: dimensionless strength of the nonlinear insert
      :param float cn: scale parameter of the nonlinear insert (m^[1/2])

   .. py:property:: particle_lost_diagnostics_backend

      Diagnostics for particles lost in apertures.
      See the ``BeamMonitor`` element for backend values.

   .. py:method:: init_grids()

      Initialize AMReX blocks/grids for domain decomposition & space charge mesh.

      This must come first, before particle beams and lattice elements are initialized.

   .. py:method:: add_particles(charge_C, distr, npart)

      Generate and add n particles to the particle container.
      Note: Set the reference particle properties (charge, mass, energy) first.

      Will also resize the geometry based on the updated particle distribution's extent and then redistribute particles in according AMReX grid boxes.

      :param float charge_C: bunch charge (C)
      :param distr: distribution function to draw from (object from :py:mod:`impactx.distribution`)
      :param int npart: number of particles to draw

   .. py:method:: particle_container()

      Access the beam particle container (:py:class:`impactx.ParticleContainer`).

   .. py:property:: lattice

      Access the elements in the accelerator lattice.
      See :py:mod:`impactx.elements` for lattice elements.

   .. py:property:: periods

      The number of periods to repeat the lattice.


   .. py:property:: abort_on_warning_threshold

      (optional) Set to "low", "medium" or "high".
      Cause the code to abort if a warning is raised that exceeds the warning threshold.

   .. py:property:: abort_on_unused_inputs

      Set to ``1`` to cause the simulation to fail *after* its completion if there were unused parameters. (default: ``0`` for false)
      It is mainly intended for continuous integration and automated testing to check that all tests and inputs are adapted to API changes.

   .. py:property:: always_warn_immediately

      If set to ``1``, ImpactX immediately prints every warning message as soon as it is generated. (default: ``0`` for false)
      It is mainly intended for debug purposes, in case a simulation crashes before a global warning report can be printed.

   .. py:method:: evolve()

      Run the main simulation loop for a number of steps.

   .. py:method:: resize_mesh()

      Resize the mesh :py:attr:`~domain` based on the :py:attr:`~dynamic_size` and related parameters.


.. py:class:: impactx.Config

      Configuration information on ImpactX that were set at compile-time.

   .. py:property:: have_mpi

      Indicates multi-process/multi-node support via the `message-passing interface (MPI) <https://www.mpi-forum.org>`__.
      Possible values: ``True``/``False``

      .. note::

         Particle beam particles are not yet dynamically load balanced.
         Please see the progress in `issue 198 <https://github.com/ECP-WarpX/impactx/issues/198>`__.

   .. py:property:: have_gpu

      Indicates GPU support.
      Possible values: ``True``/``False``

   .. py:property:: gpu_backend

      Indicates the available GPU support.
      Possible values: ``None``, ``"CUDA"`` (for Nvidia GPUs), ``"HIP"`` (for AMD GPUs) or ``"SYCL"`` (for Intel GPUs).

   .. py:property:: have_omp

      Indicates multi-threaded CPU support via `OpenMP <https://www.openmp.org>`__.
      Possible values: ``True``/``False```

      Set the environment variable ``OMP_NUM_THREADS`` to control the number of threads.

      .. warning::

         By default, OpenMP spawns as many threads as there are available virtual cores on a host.
         When MPI and OpenMP support are used at the same time, it can easily happen that one over-subscribes the available physical CPU cores.
         This will lead to a severe slow-down of the simulation.

         By setting appropriate `environment variables for OpenMP <https://www.openmp.org/spec-html/5.0/openmpch6.html>`__, ensure that the number of MPI processes (ranks) per node multiplied with the number of OpenMP threads is equal to the number of physical (or virtual) CPU cores.
         Please see our examples in the :ref:`high-performance computing (HPC) <install-hpc>` on how to run efficiently in parallel environments such as supercomputers.


Particles
---------

.. py:class:: impactx.ParticleContainer

   Beam Particles in ImpactX.

   This class stores particles, distributed over MPI ranks.

   .. py:method:: add_n_particles(x, y, t, px, py, pt, qm, bchchg)

      Add new particles to the container for fixed s.

      Note: This can only be used *after* the initialization (grids) have
            been created, meaning after the call to :py:meth:`ImpactX.init_grids`
            has been made in the ImpactX class.

      :param x: positions in x
      :param y: positions in y
      :param t: positions as time-of-flight in c*t
      :param px: momentum in x
      :param py: momentum in y
      :param pt: momentum in t
      :param qm: charge over mass in 1/eV
      :param bchchg: total charge within a bunch in C

   .. py:method:: ref_particle()

      Access the reference particle (:py:class:`impactx.RefPart`).

      :return: return a data reference to the reference particle
      :rtype: impactx.RefPart

   .. py:method:: set_ref_particle(refpart)

      Set reference particle attributes.

      :param impactx.RefPart refpart: a reference particle to copy all attributes from

   .. py:method:: reduced_beam_characteristics()

      Compute reduced beam characteristics like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters.

      :return: beam properties with string keywords
      :rtype: dict

   .. py:method:: min_and_max_positions()

      Compute the min and max of the particle position in each dimension.

      :return: x_min, y_min, z_min, x_max, y_max, z_max
      :rtype: Tuple[float, float, float, float, float, float]

   .. py:method:: mean_and_std_positions()

      Compute the mean and std of the particle position in each dimension.

      :return: x_mean, x_std, y_mean, y_std, z_mean, z_std
      :rtype: Tuple[float, float, float, float, float, float]

   .. py:method:: redistribute()

      Redistribute particles in the current mesh in x, y, z.


.. py:class:: impactx.RefPart

   This struct stores the reference particle attributes stored in :py:class:`impactx.ParticleContainer`.

   .. py:property:: s

      integrated orbit path length, in meters

   .. py:property:: x

      horizontal position x, in meters

   .. py:property:: y

      vertical position y, in meters

   .. py:property:: z

      longitudinal position y, in meters

   .. py:property:: t

      clock time * c in meters

   .. py:property:: px

      momentum in x, normalized to mass*c, :math:`p_x = \gamma \beta_x`

   .. py:property:: py

      momentum in y, normalized to mass*c, :math:`p_x = \gamma \beta_x`

   .. py:property:: pz

      momentum in z, normalized to mass*c, :math:`p_x = \gamma \beta_x`

   .. py:property:: pt

      energy, normalized by rest energy, :math:`p_t = -\gamma`

   .. py:property:: gamma

      Read-only: Get reference particle relativistic gamma, :math:`\gamma = 1/\sqrt{1-\beta^2}`

   .. py:property:: beta

      Read-only: Get reference particle relativistic beta, :math:`\beta = v/c`

   .. py:property:: beta_gamma

      Read-only: Get reference particle :math:`\beta \cdot \gamma`

   .. py:property:: qm_qeeV

      Read-only: Get reference particle charge to mass ratio (elementary charge/eV)

   .. py:method:: set_charge_qe(charge_qe)

      Write-only: Set reference particle charge in (positive) elementary charges.

   .. py:method:: set_mass_MeV(massE)

      Write-only: Set reference particle rest mass (MeV/c^2).

   .. py:method:: set_kin_energy_MeV(kin_energy_MeV)

      Write-only: Set reference particle kinetic energy (MeV)

   .. py:method:: load_file(madx_file)

      Load reference particle information from a MAD-X file.

      :param madx_file: file name to MAD-X file with a ``BEAM`` entry


Initial Beam Distributions
--------------------------

This module provides particle beam distributions that can be used to initialize particle beams in an :py:class:`impactx.ParticleContainer`.

.. py:module:: impactx.distribution
   :synopsis: Particle beam distributions in ImpactX

.. py:class:: impactx.distribution.Gaussian(sigx, sigy, sigt, sigpx, sigpy, sigpt, muxpx=0.0, muypy=0.0, mutpt=0.0)

   A 6D Gaussian distribution.

   :param sigx: for zero correlation, these are the related RMS sizes (in meters)
   :param sigy: see sigx
   :param sigt: see sigx
   :param sigpx: RMS momentum
   :param sigpy: see sigpx
   :param sigpt: see sigpx
   :param muxpx: correlation length-momentum
   :param muypy: see muxpx
   :param mutpt: see muxpx

.. py:class:: impactx.distribution.Kurth4D(sigx, sigy, sigt, sigpx, sigpy, sigpt, muxpx=0.0, muypy=0.0, mutpt=0.0)

   A 4D Kurth distribution transversely + a uniform distribution
   in t + a Gaussian distribution in pt.

.. py:class:: impactx.distribution.Kurth6D(sigx, sigy, sigt, sigpx, sigpy, sigpt, muxpx=0.0, muypy=0.0, mutpt=0.0)

   A 6D Kurth distribution.

   R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)
   C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)

.. py:class:: impactx.distribution.KVdist(sigx, sigy, sigt, sigpx, sigpy, sigpt, muxpx=0.0, muypy=0.0, mutpt=0.0)

   A K-V distribution transversely + a uniform distribution
   in t + a Gaussian distribution in pt.

.. py:class:: impactx.distribution.None

   This distribution sets all values to zero.

.. py:class:: impactx.distribution.Semigaussian(sigx, sigy, sigt, sigpx, sigpy, sigpt, muxpx=0.0, muypy=0.0, mutpt=0.0)

   A 6D Semi-Gaussian distribution (uniform in position, Gaussian in momentum).

.. py:class:: impactx.distribution.Triangle(sigx, sigy, sigt, sigpx, sigpy, sigpt, muxpx=0.0, muypy=0.0, mutpt=0.0)

   A triangle distribution for laser-plasma acceleration related applications.

   A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).
   The transverse distribution is a 4D waterbag.

.. py:class:: impactx.distribution.Waterbag(sigx, sigy, sigt, sigpx, sigpy, sigpt, muxpx=0.0, muypy=0.0, mutpt=0.0)

   A 6D Waterbag distribution.

.. py:class:: impactx.distribution.Thermal(k, kT, kT_halo, normalize, normalize_halo, halo)

   A 6D stationary thermal or bithermal distribution.


Lattice Elements
----------------

This module provides elements for the accelerator lattice.

.. py:module:: impactx.elements
   :synopsis: Accelerator lattice elements in ImpactX

.. py:class:: impactx.elements.KnownElementsList

   An iterable, ``list``-like type of elements.

   .. py:method:: clear()

      Clear the list to become empty.

   .. py:method:: extend(list)

      Add a list of elements to the list.

   .. py:method:: append(element)

      Add a single element to the list.

   .. py:method:: load_file(madx_file, nslice=1)

      Load and append an accelerator lattice description from a MAD-X file.

      :param madx_file: file name to MAD-X file with beamline elements
      :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.CFbend(ds, rc, k, dx=0, dy=0, rotation=0, nslice=1)

   A combined function bending magnet.  This is an ideal Sbend with a normal quadrupole field component.

   :param ds: Segment length in m.
   :param rc: Radius of curvature in m.
   :param k:  Quadrupole strength in m^(-2) (MADX convention)
              = (gradient in T/m) / (rigidity in T-m)
              k > 0 horizontal focusing
              k < 0 horizontal defocusing
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.ConstF(ds, kx, ky, kt, dx=0, dy=0, rotation=0, nslice=1)

   A linear Constant Focusing element.

   :param ds: Segment length in m.
   :param kx: Focusing strength for x in 1/m.
   :param ky: Focusing strength for y in 1/m.
   :param kt: Focusing strength for t in 1/m.
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

   .. py:property:: kx

      focusing x strength in 1/m

   .. py:property:: ky

      focusing y strength in 1/m

   .. py:property:: kt

      focusing t strength in 1/m

.. py:class:: impactx.elements.DipEdge(psi, rc, g, K2, dx=0, dy=0, rotation=0)

   Edge focusing associated with bend entry or exit

   This model assumes a first-order effect of nonzero gap.
   Here we use the linear fringe field map, given to first order in g/rc (gap / radius of curvature).

   References:

   * K. L. Brown, SLAC Report No. 75 (1982).
   * K. Hwang and S. Y. Lee, PRAB 18, 122401 (2015).

   :param psi: Pole face angle in rad
   :param rc: Radius of curvature in m
   :param g: Gap parameter in m
   :param K2: Fringe field integral (unitless)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

.. py:class:: impactx.elements.Drift(ds, dx=0, dy=0, rotation=0, nslice=1)

   A drift.

   :param ds: Segment length in m
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.ChrDrift(ds, dx=0, dy=0, rotation=0, nslice=1)

   A drift with chromatic effects included.  The Hamiltonian is expanded
   through second order in the transverse variables (x,px,y,py), with the exact pt
   dependence retained.

   :param ds: Segment length in m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.ExactDrift(ds, dx=0, dy=0, rotation=0, nslice=1)

   A drift using the exact nonlinear transfer map.

   :param ds: Segment length in m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.Kicker(xkick, ykick, units, dx=0, dy=0, rotation=0)

   A thin transverse kicker.

   :param xkick: horizontal kick strength (dimensionless OR T-m)
   :param ykick: vertical kick strength (dimensionless OR T-m)
   :param units: specification of units (``"dimensionless"`` in units of the magnetic rigidity of the reference particle or ``"T-m"``)

.. py:class:: impactx.elements.Multipole(multipole, K_normal, K_skew, dx=0, dy=0, rotation=0)

   A general thin multipole element.

   :param multipole: index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
   :param K_normal: Integrated normal multipole coefficient (1/meter^m)
   :param K_skew: Integrated skew multipole coefficient (1/meter^m)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

.. py::class:: impactx.elements.None

   This element does nothing.

.. py:class:: impactx.elements.NonlinearLens(knll, cnll, dx=0, dy=0, rotation=0)

   Single short segment of the nonlinear magnetic insert element.

   A thin lens associated with a single short segment of the
   nonlinear magnetic insert described by V. Danilov and
   S. Nagaitsev, PRSTAB 13, 084002 (2010), Sect. V.A.  This
   element appears in MAD-X as type ``NLLENS``.

   :param knll: integrated strength of the nonlinear lens (m)
   :param cnll: distance of singularities from the origin (m)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

.. py:class:: impactx.elements.BeamMonitor(name, backend="default", encoding="g")

   A beam monitor, writing all beam particles at fixed ``s`` to openPMD files.

   If the same element ``name`` is used multiple times, then an output series is created with multiple outputs.

   The `I/O backend <https://openpmd-api.readthedocs.io/en/latest/backends/overview.html>`_ for `openPMD <https://www.openPMD.org>`_ data dumps.
   ``bp`` is the `ADIOS2 I/O library <https://csmd.ornl.gov/adios>`_, ``h5`` is the `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_, and ``json`` is a `simple text format <https://en.wikipedia.org/wiki/JSON>`_.
   ``json`` only works with serial/single-rank jobs.
   By default, the first available backend in the order given above is taken.

   openPMD `iteration encoding <https://openpmd-api.readthedocs.io/en/0.14.0/usage/concepts.html#iteration-and-series>`__ determines if multiple files are created for individual output steps or not.
   Variable based is an `experimental feature with ADIOS2 <https://openpmd-api.readthedocs.io/en/0.14.0/backends/adios2.html#experimental-new-adios2-schema>`__.

   :param name: name of the series
   :param backend: I/O backend, e.g., ``bp``, ``h5``, ``json``
   :param encoding: openPMD iteration encoding: (v)ariable based, (f)ile based, (g)roup based (default)

.. py:class:: impactx.elements.Programmable

   A programmable beam optics element.

   This element can be programmed to receive callback hooks into Python functions.

   .. py:property:: beam_particles

      This is a function hook for pushing all beam particles.
      This accepts a function or lambda with the following arguments:

      .. py:method:: user_defined_function(pti: ImpactXParIter, refpart: RefPart)

         This function is called repeatedly for all particle tiles or boxes in the beam particle container.
         Particles can be pushed and are relative to the reference particle

   .. py:property:: ref_particle

      This is a function hook for pushing the reference particle.
      This accepts a function or lambda with the following argument:

      .. py:method:: another_user_defined_function(refpart: RefPart)

         This function is called for the reference particle as it passes through the element.
         The reference particle is updated *before* the beam particles are pushed.

.. py:class:: impactx.elements.Quad(ds, k, dx=0, dy=0, rotation=0, nslice=1)

   A Quadrupole magnet.

   :param ds: Segment length in m.
   :param k:  Quadrupole strength in m^(-2) (MADX convention)
              = (gradient in T/m) / (rigidity in T-m)
              k > 0 horizontal focusing
              k < 0 horizontal defocusing
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.ChrQuad(ds, k, units, dx=0, dy=0, rotation=0, nslice=1)

   A Quadrupole magnet, with chromatic effects included.  The Hamiltonian is expanded
   through second order in the transverse variables (x,px,y,py), with the exact pt
   dependence retained.

   :param ds: Segment length in m.
   :param k:  Quadrupole strength in m^(-2) (MADX convention, if units = 0)
              = (gradient in T/m) / (rigidity in T-m)
          OR  Quadrupole strength in T/m (MaryLie convention, if units = 1)
              k > 0 horizontal focusing
              k < 0 horizontal defocusing
   :param units: specification of units for quadrupole field strength
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

   .. py:property:: k

      quadrupole strength in 1/m^2 (or T/m)

   .. py:property:: units

      unit specification for quad strength

.. py:class:: impactx.elements.ChrAcc(ds, ez, bz, dx=0, dy=0, rotation=0, nslice=1)

   Acceleration in a uniform field Ez, with a uniform solenoidal field Bz.

   The Hamiltonian is expanded through second order in the
   transverse variables (x,px,y,py), with the exact pt dependence retained.

   :param ds: Segment length in m
   :param ez: electric field strength in m^(-1)
              = (charge * electric field Ez in V/m) / (m*c^2)
   :param bz: magnetic field strength in m^(-1)
              = (charge * magnetic field Bz in T) / (m*c)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

   .. py:property:: ez

      electric field strength in 1/m

   .. py:property:: bz

      magnetic field strength in 1/m

.. py:class:: impactx.elements.RFCavity(ds, escale, freq, phase, dx=0, dy=0, rotation=0, mapsteps=1, nslice=1)

   A radiofrequency cavity.

   :param ds: Segment length in m.
   :param escale: scaling factor for on-axis RF electric field in 1/m
                  = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV)
   :param freq: RF frequency in Hz
   :param phase: RF driven phase in degrees
   :param cos_coefficients: array of ``float`` cosine coefficients in Fourier expansion of on-axis electric field Ez (optional); default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__

   :param cos_coefficients: array of ``float`` sine coefficients in Fourier expansion of on-axis electric field Ez (optional); default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param mapsteps: number of integration steps per slice used for map and reference particle push in applied fields
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.Sbend(ds, rc, dx=0, dy=0, rotation=0, nslice=1)

   An ideal sector bend.

   :param ds: Segment length in m.
   :param rc: Radius of curvature in m.
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.ExactSbend(ds, phi, B, dx=0, dy=0, rotation=0, nslice=1)

   An ideal sector bend using the exact nonlinear map.  The model consists of a uniform bending field B_y with a hard edge.  Pole faces are
   normal to the entry and exit velocity of the reference particle.

   References:

   * D. L. Bruhwiler et al, in Proc. of EPAC 98, pp. 1171-1173 (1998).
   * E. Forest et al, Part. Accel. 45, pp. 65-94 (1994).

   :param ds: Segment length in m.
   :param phi: Bend angle in degrees.
   :param B: Magnetic field in Tesla; when B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad),   corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B.
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.Buncher(V, k, dx=0, dy=0, rotation=0)

   A short RF cavity element at zero crossing for bunching (MaryLie model).

   :param V: Normalized RF voltage drop V = Emax*L/(c*Brho)
   :param k: Wavenumber of RF in 1/m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

.. py:class:: impactx.elements.ShortRF(V, freq, phase, dx=0, dy=0, rotation=0)

   A short RF cavity element (MAD-X model).

   :param V: Normalized RF voltage V = maximum energy gain/(m*c^2)
   :param freq: RF frequency in Hz
   :param phase: RF synchronous phase in degrees (phase = 0 corresponds to maximum energy gain, phase = -90 corresponds go zero energy gain for bunching)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

.. py:class:: impactx.elements.ChrUniformAcc(ds, k, dx=0, dy=0, rotation=0, nslice=1)

   A region of constant Ez and Bz for uniform acceleration, with chromatic effects included.
   The Hamiltonian is expanded through second order in the transverse variables (x,px,y,py),
   with the exact pt dependence retained.

   :param ds: Segment length in m.
   :param ez: Electric field strength in m^(-1)
              = (particle charge in C * field Ez in V/m) / (particle mass in kg * (speed of light in m/s)^2)
   :param bz: Magnetic field strength in m^(-1)
              = (particle charge in C * field Bz in T) / (particle mass in kg * speed of light in m/s)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.SoftSolenoid(ds, bscale, cos_coefficients, sin_coefficients, dx=0, dy=0, rotation=0, mapsteps=1, nslice=1)

   A soft-edge solenoid.

   :param ds: Segment length in m.
   :param bscale: Scaling factor for on-axis magnetic field Bz in inverse meters
   :param cos_coefficients: array of ``float`` cosine coefficients in Fourier expansion of on-axis magnetic field Bz
            (optional); default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__
   :param sin_coefficients: array of ``float`` sine coefficients in Fourier expansion of on-axis magnetic field Bz
            (optional); default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param mapsteps: number of integration steps per slice used for map and reference particle push in applied fields
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.Sol(ds, ks, dx=0, dy=0, rotation=0, nslice=1)

   An ideal hard-edge Solenoid magnet.

   :param ds: Segment length in m.
   :param ks: Solenoid strength in m^(-1) (MADX convention) in (magnetic field Bz in T) / (rigidity in T-m)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.PRot(phi_in, phi_out)

   Exact map for a pole-face rotation in the x-z plane.

   :param phi_in: angle of the reference particle with respect to the longitudinal (z) axis in the original frame in degrees
   :param phi_out: angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in degrees
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

.. py:class:: impactx.elements.Aperture(xmax, ymax, shape="rectangular", dx=0, dy=0, rotation=0)

   A thin collimator element, applying a transverse aperture boundary.

   :param xmax: maximum allowed value of the horizontal coordinate (meter)
   :param ymax: maximum allowed value of the vertical coordinate (meter)
   :param shape: aperture boundary shape: ``"rectangular"`` (default) or ``"elliptical"``
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

   .. py:property:: shape

      aperture type (rectangular, elliptical)

   .. py:property:: xmax

      maximum horizontal coordinate

   .. py:property:: ymax

      maximum vertical coordinate

.. py:class:: impactx.elements.SoftQuadrupole(ds, gscale, cos_coefficients, sin_coefficients, dx=0, dy=0, rotation=0, mapsteps=1, nslice=1)

   A soft-edge quadrupole.

   :param ds: Segment length in m.
   :param gscale: Scaling factor for on-axis field gradient in inverse meters
   :param cos_coefficients: array of ``float`` cosine coefficients in Fourier expansion of on-axis field gradient
            (optional); default is a tanh fringe field model based on `<http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__
   :param sin_coefficients: array of ``float`` sine coefficients in Fourier expansion of on-axis field gradient
            (optional); default is a tanh fringe field model based on `<http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param mapsteps: number of integration steps per slice used for map and reference particle push in applied fields
   :param nslice: number of slices used for the application of space charge

.. py:class:: impactx.elements.ThinDipole(theta, rc, dx=0, dy=0, rotation=0)

   A general thin dipole element.

   :param theta: Bend angle (degrees)
   :param rc: Effective curvature radius (meters)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]

   Reference:

   * G. Ripken and F. Schmidt, Thin-Lens Formalism for Tracking, CERN/SL/95-12 (AP), 1995.


Coordinate Transformation
-------------------------

.. py:class:: impactx.TransformationDirection

   Enumerated type indicating whether to transform to fixed :math:`s` or fixed :math:`t` coordinate system when applying ``impactx.coordinate_transformation``.

   :param to_fixed_t:
   :param to_fixed_s:

.. py:function:: impactx.coordinate_transformation(pc, direction)

   Function to transform the coordinates of the particles in a particle container either to fixed :math:`t` or to fixed :math:`s`.

   :param pc: ``impactx.particle_container`` whose particle coordinates are to be transformed.
   :param direction: enumerated type ``impactx.TransformationDirection``, indicates whether to transform to fixed :math:`s` or fixed :math:`t`.
