.. _usage-python:

Parameters: Python
==================

This documents on how to use ImpactX as a Python script (``python3 run_script.py``).

.. tip::

   If you enjoy AI/LLM/agentic workflows, see our :ref:`AI (LLM)-Assisted Input File Design <ai_input_design>` section, too.

Collective Effects & Overall Simulation Parameters
--------------------------------------------------

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

   .. py:property:: max_grid_size

      The maximum number of grid points per box (grid) along each direction, as a list with one entry per mesh-refinement level.
      The mesh is decomposed into boxes of at most this size.
      For best performance, aim for exactly one box per parallel (MPI) process, which is also one box per GPU on accelerated machines; ImpactX warns at the start of a run if the number of boxes does not match the number of processes.
      A direction that is also set through :py:attr:`~max_grid_size_x`, :py:attr:`~max_grid_size_y`, :py:attr:`~max_grid_size_z` or :py:attr:`~blocking_factor_x` keeps the per-direction value, no matter in which order the properties are assigned.
      Use the per-direction properties in that case.

   .. py:property:: max_grid_size_x

      The maximum number of grid points per box along ``x``, as a list with one entry per mesh-refinement level, overwriting :py:attr:`~max_grid_size`.
      This property is also set when assigning :py:attr:`~blocking_factor_x`.
      The per-direction value always wins over :py:attr:`~max_grid_size`, independent of the order of assignment, so assign this property to change ``x``.
      Analogous properties exist as ``max_grid_size_y`` and ``max_grid_size_z``.

   .. py:property:: blocking_factor_x

      The AMReX blocking factor for the ``x`` direction, as a list with one entry per mesh-refinement level.
      Box sizes are forced to be a multiple of this value, which must be a power of two.
      Setting this property also sets :py:attr:`~max_grid_size_x` to the same value, which preserves one box per process in the default decomposition; assign :py:attr:`~max_grid_size_x` afterwards to overwrite.
      Assigning the all-direction :py:attr:`~max_grid_size` does not overwrite it, since per-direction values always win.
      Analogous properties exist as ``blocking_factor_y`` and ``blocking_factor_z``.

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

      This is a list with :pp:param:`amr.max_level` + 1 entries.

      By default, we dynamically extract the minimum and maximum of the particle positions in the beam.
      The field mesh spans, per direction, multiple times the maximum physical extent of beam particles, as given by this factor.
      The beam minimum and maximum extent are symmetrically padded by the mesh.
      For instance, ``1.2`` means the mesh will span 10% above and 10% below the beam;
      ``1.0`` means the beam is exactly covered with the mesh.

      Default: ``[3.0 1.0 1.0 ...]``.
      When set, turns ``dynamic_size`` to ``True``.

   .. py:property:: dynamic_size

      Use dynamic (``True``) resizing of the field mesh or static sizing (``False``).

   .. py:property:: strang_split

      Compose the collective effect kicks and the element transport in a second-order,
      time-symmetric Strang split (``True``, default) or to first order (``False``).
      See :pp:param:`algo.strang_split`.

   .. py:property:: space_charge

      The physical model of space charge used.

      Options:

      * ``False`` (default): space charge effects are not calculated.

      * ``"2D"``: Space charge forces are computed in the plane ``(x,y)`` transverse to the reference particle velocity, assuming the beam is long and unbunched.

      * ``"2p5D"``: Space charge forces are computed in the plane ``(x,y)`` transverse to the reference particle velocity, while the transverse space charge kicks are weighted by the
        longitudinal line density determined by charge deposition (2.5D model).  Longitudinal space charge kicks are determined by the derivative of the line charge density.

      * ``"3D"``: Space charge forces are computed in three dimensions, assuming the beam is bunched.

        When running in envelope mode (when :pp:param:`algo.track = envelope`), this model currently assumes that ``<xy> = <yt> = <tx> = 0``.

      * ``"Gauss3D"``: Calculate 3D space charge forces as if the beam was a Gaussian distribution.

      * ``"Gauss2p5D"``: Calculate 2.5D space charge forces as if the beam was a transverse Gaussian distribution.

        These models are supported only in particle tracking mode (when :pp:param:`algo.track = particles`).
        Ref.: J. Qiang, "Two-and-a-half dimensional symplectic space-charge solver", LBNL Report Number: LBNL-2001674 (2025).
        (This reference describes both 3D and 2.5D models.)

   .. py:property:: space_charge_gauss_nint

      Number of steps for computing the integrals (default: ``101``).

   .. py:property:: space_charge_gauss_taylor_delta

      Initial integral region to avoid integrand divergence at 0 (default: ``0.01``).

   .. py:property:: space_charge_gauss_charge_z_bins

      Number of bins for longitudinal charge density deposition (default: ``129``).  Used by the Gauss2p5D space charge model.

   .. py:property:: space_charge_gauss_long_scale

      Longitudinal space charge scale for the Gauss2p5D space charge model.
      This is an approximation that only influences the longitudinal momentum (``pt``) kick.
      If not set, it defaults to :math:`1.103 \cdot \gamma \cdot \sigma_z`, estimated in-situ from the
      current reduced beam characteristics (with :math:`\sigma_z` the RMS bunch length), which is a
      typical value when comparing to a 3D model.  This value minimizes the L2-norm of the error in
      the on-axis field Ez in the case of a long 3D Gaussian bunch.
      In a perfectly conducting pipe with a long bunch, where the bunch length in the
      rest frame is significantly larger than the pipe radius, it is recommended to set
      the longitudinal scale length to the pipe radius in the input file.

   .. py:property:: space_charge_num_longitudinal_bins

      Number of bins for longitudinal charge density deposition (default: ``100``).  Used by the 2p5D space charge model.

   .. py:property:: space_charge_apply_longitudinal_kick

      Enable or disable the longitudinal space charge kick (default: ``True``).

   .. py:property:: poisson_solver

      The numerical solver to solve the Poisson equation when calculating space charge effects.
      Either ``"fft"`` (default) or ``"multigrid"``.

      Currently, the multigrid solver supports only 3D space charge.  The fft solver supports 2D, 2.5D or 3D space charge.

      * ``fft``: Poisson's equation is solved using an Integrated Green Function method (which requires FFT calculations).
        See these references for more details `Qiang et al. (2006) <https://doi.org/10.1103/PhysRevSTAB.9.044204>`__ (+ `Erratum <https://doi.org/10.1103/PhysRevSTAB.10.129901>`__).
        This requires the compilation flag ``-DImpactX_FFT=ON``.
        If mesh refinement (MR) is enabled, this FFT solver is used only on the coarsest level and a multi-grid solver is used on refined levels.
        The boundary conditions are assumed to be open.

      * ``multigrid``: Poisson's equation is solved using an iterative multigrid (MLMG) solver.
        See the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#>`__ for details of the MLMG solver.

   .. py:property:: mlmg_relative_tolerance

      Default: ``1.e-7`` (DP) / ``1.e-4`` (SP)

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

   .. py:property:: csr

      Enable (``True``) or disable (``False``) Coherent Synchrotron Radiation (CSR) calculations (default: ``False``).

      Whether to calculate Coherent Synchrotron Radiation (CSR) effects (default: disabled).

      Currently, this is the 1D ultrarelativistic steady-state wakefield model (eq. 19 of
      `E. L. Saldin et al, NIMA 398, p. 373-394 (1997), DOI:10.1016/S0168-9002(97)00822-X <https://doi.org/10.1016/S0168-9002(97)00822-X>`__).

      .. note::

         CSR effects are only calculated for lattice elements that include bending, such as ``Sbend``, ``ExactSbend`` and ``CFbend``.

         CSR effects require the compilation flag ``-DImpactX_FFT=ON``.

   .. py:property:: csr_bins

      The number of bins along the longitudinal direction used for the CSR calculations (default: ``150``).

   .. py:property:: isr

      Enable (``True``) or disable (``False``) Incoherent Synchrotron Radiation (ISR) calculations (default: ``False``).

      Whether to calculate Incoherent Synchrotron Radiation (ISR) effects (default: disabled).

      ISR effects are included in the simulation for bend lattice elements such as ``Sbend`` and ``CFbend``, and are especially important for electron or positron bunches at high energy.
      The effects of ISR include radiation reaction due to the stochastic emission of synchrotron radiation, resulting in mean energy loss and quantum excitation of the bunch.
      The model is based on:

      `F. Niel et al., Phys. Rev. E 97, 043209 (2018), DOI:10.1103/PhysRevE.97.043209 <https://doi.org/10.1103/PhysRevE.97.043209>`__

      However, a Taylor expansion is used to evaluate the dependence on the quantum parameter :math:`\chi`.  When :pp:param:`algo.isr_order = 1`, the model is equivalent to that described in:

      `J. M. Jowett, "Introductory Statistical Mechanics for Electron Storage Rings", AIP Conf. Proc. 153, 864-970 (1987), DOI:10.1063/1.36374 <https://doi.org/10.1063/1.36374>`__

      .. note::

         ISR effects are only calculated for lattice elements that include bending, such as ``Sbend``, ``ExactSbend`` and ``CFbend``.

   .. py:property:: isr_order

      The number of terms retained in the Taylor series for the functions :math:`g(\chi)` and :math:`h(\chi)` appearing in Niel et al, equations (25) and (41) describing quantum effects.

   .. py:property:: isr_on_ref_part

      Flag specifying whether ISR is to be applied to the reference particle.  When ``sim.isr_on_ref_part = False``, the reference particle does not lose energy due to radiation, and the
      mean energy of the beam particles will decrease.  This option is natural if the lattice optics, magnet settings, etc. are chosen without accounting for radiative energy loss.
      When ``sim.isr_on_ref_part = True``, the reference particle does lose energy due to radiation, and little centroid evolution is expected in the beam particles.  This option is natural if the lattice optics, magnet settings, etc. are chosen to account for radiative energy loss.

   .. py:property:: particle_bc

      The application of a non-trivial boundary condition for particles is currently supported only in the longitudinal direction (the local direction of motion as defined by the reference
      particle).  This parameter sets the type of particle boundary condition to be applied to the longitudinal coordinate, based on the value of parameter ``bucket_length``.

      Options:

      * ``"open"`` (default): no action is applied at the boundary.

      * ``"periodic"``: each particle's longitudinal coordinate is treated modulo the value ``bucket_length`` (phase wrapping).

      * ``"absorbing"``: a particle whose longitudinal coordinate falls outside the boundary is declared lost.

      * ``"reflecting"``: a particle whose longitudinal coordinate crosses the boundary is reflected about the boundary, with reversed longitudinal momentum.

        .. note::

           The implementation works through linear order in the phase space variables.
           If you have the need for a more precise implementation of reflecting boundaries, please `open an issue <https://github.com/BLAST-ImpactX/impactx/issues/new>`__.

   .. py:property:: spin

      Enable (``True``) or disable (``False``) particle spin tracking (default: ``False``).

      Whether to track particle spin.
      Currently, the implementation of spin tracking is a work in progress, and this feature is not yet supported.

      Spin tracking uses the gyromagnetic anomaly of the reference particle.
      Set it via ``RefPart.set_species()`` or ``RefPart.set_gyromagnetic_anomaly()``.

   .. py:property:: diagnostics

      Enable (``True``) or disable (``False``) diagnostics generally (default: ``True``).
      Disabling this is mostly used for benchmarking.

   .. py:property:: slice_step_diagnostics

      Enable (``True``) or disable (``False``) diagnostics every slice step in elements  (default: ``True``).

      By default, diagnostics is performed at the beginning and end of the simulation.
      Enabling this flag will write diagnostics every step and slice step.

   .. py:property:: diag_file_prefix

      Root directory for diagnostic output (default: folder named ``"diags"``).

      Set to ``""`` or ``"."`` to write diagnostics in the current working directory.

      If a directory at ``diag_file_prefix`` already exists when a simulation starts,
      ImpactX renames it to ``<diag_file_prefix>.old.<suffix>`` to preserve prior results.
      This is skipped when ``diag_file_prefix`` resolves to the current working directory,
      the root directory, or an ancestor of the current working directory; in those cases
      new output is written alongside existing files.

   .. py:property:: diag_file_min_digits

      The minimum number of digits (default: ``6``) used for the step
      number appended to the diagnostic file names.

   .. py:property:: particle_lost_diagnostics_backend

      Diagnostics for particles lost in apertures.
      See the ``BeamMonitor`` element for backend values.

   .. py:property:: eigenemittances

      Enable (``True``) or disable (``False``) output of eigenemittances at every slice step in elements  (default: ``False``).

      If this flag is enabled, the 3 eigenemittances of the 6D beam distribution are computed and written as diagnostics.
      This flag is disabled by default to reduce computational cost.

   .. py:method:: init_grids()

      Initialize AMReX blocks/grids for domain decomposition & space charge mesh.

      This must come first, before particle beams and lattice elements are initialized.

   .. py:method:: add_particles(charge_C, distr, npart, spinv=None)

      Particle tracking mode: Generate and add n particles to the particle container.
      Note: Set the reference particle properties (charge, mass, energy) first.

      Will also resize the geometry based on the updated particle distribution's extent and then redistribute particles in according AMReX grid boxes.

      :param float charge_C: bunch charge (C)
      :param distr: distribution function to draw from (object from :py:mod:`impactx.distribution`)
      :param int npart: number of particles to draw
      :param SpinvMF spinv: optional spin distribution

   .. py:method:: init_envelope(ref, distr, intensity=None)

      Envelope tracking mode:
      Create a 6x6 covariance matrix from a distribution and then initialize the simulation for envelope tracking relative to a reference particle.

      :param ref: the reference particle (object from :py:class:`impactx.RefPart`)
      :param distr: distribution function (object from :py:mod:`impactx.distribution`)
      :param float intensity: the beam intensity, given as bunch charge (C) for 3D or beam current (A) for 2D space charge

   .. py:method:: particle_container()

      Access the beam particle container (:py:class:`impactx.ParticleContainer`).
      Deprecated: use :py:attr:`beam`.

   .. py:property:: beam

      Access the beam particle container (:py:class:`impactx.ParticleContainer`).

      For example, ``sim.beam.to_df(local=True)`` returns the local particles as a pandas DataFrame.
      See :ref:`usage-howto-python-particle-data` for further data-access recipes.

   .. py:property:: envelope

      Access the beam envelope (:py:class:`impactx.Envelope`: a 6x6 covariance matrix
      and beam intensity) used for envelope tracking.
      Only available after :py:meth:`init_envelope`.

   .. py:method:: rho(lev)

      Return the charge density :math:`\rho` on mesh-refinement level ``lev`` as a pyAMReX ``MultiFab``.
      Populated when space-charge is enabled (see :py:attr:`space_charge`).

      :param int lev: mesh-refinement level.

      .. note::

         Field data are populated by the space-charge solve at each slice step.
         They are most reliably read from inside a :py:attr:`hook` ``"after_element"`` callback
         (the element's last slice has just solved) or right after :py:meth:`track_particles` returns.
         See :ref:`usage-howto-python-field-data` and the
         `pyAMReX MultiFab guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#fields>`__
         for indexing details.

   .. py:method:: phi(lev)

      Return the scalar potential :math:`\phi` on mesh-refinement level ``lev`` as a pyAMReX ``MultiFab``.
      Populated when space-charge is enabled.
      See lifetime caveats and indexing in :ref:`usage-howto-python-field-data`.

      :param int lev: mesh-refinement level.

   .. py:method:: space_charge_field(lev, comp)

      Return one Cartesian component of the space-charge force as a pyAMReX ``MultiFab``.

      :param int lev: mesh-refinement level.
      :param str comp: field component to access (``"x"``, ``"y"``, or ``"z"``).

      See lifetime caveats and indexing in :ref:`usage-howto-python-field-data`.

   .. py:method:: Geom(lev)

      Return the AMReX ``Geometry`` object for mesh-refinement level ``lev``.
      Useful to query the physical domain extent and cell sizes for field data.

   .. py:method:: boxArray(lev)

      Return the AMReX ``BoxArray`` for mesh-refinement level ``lev``.

   .. py:method:: DistributionMap(lev)

      Return the AMReX ``DistributionMapping`` for mesh-refinement level ``lev``.

   .. py:property:: lattice

      Access the elements in the accelerator lattice.
      See :py:mod:`impactx.elements` for lattice elements.

   .. py:property:: periods

      The number of periods to repeat the lattice.

   .. py:property:: omp_threads

      Controls the number of OpenMP threads to use (ImpactX default: "nosmt").

      See the detailed `AMReX docs <https://amrex-codes.github.io/amrex/docs_html/InputsComputeBackends.html>`__ for details in the accepted values.

   .. py:property:: abort_on_warning_threshold

      (optional) Set to "low", "medium" or "high".
      Cause the code to abort if a warning is raised that exceeds the warning threshold.

   .. py:property:: abort_on_unused_inputs

      Set to ``1`` to cause the simulation to fail *after* its completion if there were unused parameters. (default: ``0`` for false)
      It is mainly intended for continuous integration and automated testing to check that all tests and inputs are adapted to API changes.

   .. py:property:: always_warn_immediately

      If set to ``1``, ImpactX immediately prints every warning message as soon as it is generated. (default: ``0`` for false)
      It is mainly intended for debug purposes, in case a simulation crashes before a global warning report can be printed.

   .. py:property:: verbose

      Controls how much information is printed to the terminal, when running ImpactX.
      ``0`` for silent, higher is more verbose. Default is ``1``.

   .. py:property:: tiny_profiler

      This parameter can be used to disable tiny profiling including CArena memory profiling at runtime.
      Default is ``True``.

   .. py:property:: memory_profiler

      This parameter can be used to disable tiny profiler's memory arena profiling at runtime.
      If ```tiny_profiler`` is ``False``, this parameter has no effects.
      Default is ``True``.

   .. py:property:: tiny_profiler_file

      If this parameter is empty (default), the output of tiny profiling is dumped on the default out stream of AMReX.
      If it's not empty, it specifies the file name for the output.
      Note that ``"/dev/null"`` is a special name that mean no output.

   .. py:method:: evolve()

      Run the main simulation loop (deprecated, use ``track_particles``)

   .. py:method:: track_particles()

      Run the particle tracking simulation loop.

   .. py:method:: track_envelope()

      Run the envelope tracking simulation loop.

      .. note::

         Our current envelope tracking implements ideal transfer maps, assuming always zero misalignments for translations.
         Element rotations are handled.
         Support for translations errors, non-zero envelope means, and feed-down effects in envelope tracking is in development.
         Until then, translations errors set on elements are silently ignored.

   .. py:method:: track_reference(ref)

      Run the reference orbit tracking simulation loop.

      :param ref: the reference particle (object from :py:class:`impactx.RefPart`)

   .. py:property:: hook

      User-defined function hooks that are called, e.g, during tracking.
      Supported hook locations names are:

      * ``"before_period"``: before each period (e.g., turn or channel period)
      * ``"after_period"``: after each period (e.g., turn or channel period)
      * ``"before_element"``: before each element is entered
      * ``"after_element"``: after each element is exited
      * ``"before_slice"``: before each element slice

      Example: Function hook that can be called before each turn (sim):

      .. code-block:: python3

         def hook_before_period(sim):
             beam = sim.beam
             turn = sim.tracking_period
             # Example: you could now manipulate elements in sim.lattice
             #          for the next turn.

         sim.hook["before_period"] = hook_before_period

      Full example: :ref:`Acceleration by RF Cavities <examples-rfcavity-ref-part-hook>`.

      See :ref:`usage-howto-python-extend` for more callback recipes.

   .. py:property:: tracking_step

      For tracking hooks/callbacks, a global step of the simulation.

      A state of internal simulation steps, increments also for space charge slice steps in elements.
      We start in "step 0" (initial state).

   .. py:property:: tracking_period

      For tracking hooks/callbacks, the period in the lattice (e.g., turn or channel period).

   .. py:property:: tracking_element

      For tracking hooks/callbacks, the current lattice element.

   .. py:method:: resize_mesh()

      Resize the mesh :py:attr:`~domain` based on the :py:attr:`~dynamic_size` and related parameters.


.. py:class:: impactx.Envelope

   The beam envelope used for envelope (covariance) tracking: a 6x6 covariance
   matrix relative to a reference particle, plus an optional beam intensity.
   Created from a distribution by :py:meth:`impactx.ImpactX.init_envelope` and
   accessed (during or after :py:meth:`impactx.ImpactX.track_envelope`) via
   :py:attr:`impactx.ImpactX.envelope`.

   .. py:property:: envelope

      The 6x6 beam covariance matrix (a ``Map6x6``).

   .. py:property:: beam_intensity

      The beam intensity: bunch charge (C) for 3D or beam current (A) for 2D space charge.

   .. py:method:: beam_moments(ref)

      Calculate beam moments (position and momentum moments, emittances, Twiss
      functions, dispersion, ...) from this envelope's covariance matrix and a
      reference particle. The envelope counterpart of
      :py:meth:`impactx.ParticleContainer.beam_moments`.

      :param ref: the reference particle (object from :py:class:`impactx.RefPart`)
      :return: a dictionary of beam moments
      :rtype: dict[str, float]


.. py:class:: impactx.Config

      Configuration information on ImpactX that were set at compile-time.

   .. py:property:: have_mpi

      Indicates multi-process/multi-node support via the `message-passing interface (MPI) <https://www.mpi-forum.org>`__.
      Possible values: ``True``/``False``

      .. note::

         Particle beam particles are not yet dynamically load balanced.
         Please see the progress in `issue 198 <https://github.com/BLAST-ImpactX/impactx/issues/198>`__.

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

         By default, OpenMP spawns as many threads as there are available physical CPU cores on a host.
         When MPI and OpenMP support are used at the same time, it can easily happen that one over-subscribes the available physical CPU cores.
         This will lead to a severe slow-down of the simulation.

         By setting appropriate `environment variables for OpenMP <https://www.openmp.org/spec-html/5.0/openmpch6.html>`__, ensure that the number of MPI processes (ranks) per node multiplied with the number of OpenMP threads is equal to the number of physical (or virtual) CPU cores.
         Please see our examples in the :ref:`high-performance computing (HPC) <install-hpc>` on how to run efficiently in parallel environments such as supercomputers.


Particles
---------

:py:class:`impactx.ParticleContainer` derives from a pyAMReX particle container.
Many bulk-access methods used to read or modify the beam in-memory (such as :py:meth:`~impactx.ParticleContainer.to_df` for a pandas DataFrame, or iterating over particle tiles) come from there.
For a full reference of the inherited compute API, see the `pyAMReX particles guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__.

For step-by-step recipes on how to access particle data live during a simulation, see :ref:`usage-howto-python-particle-data`.

.. py:class:: impactx.ParticleContainer

   Beam Particles in ImpactX.

   This class stores particles, distributed over MPI ranks.

   .. py:method:: clear(keep_mass=False, keep_charge=False)

      Empty the container and reset the reference particle.

      :param keep_mass: do not reset the reference particle mass
      :param keep_charge: do not reset the reference particle charge

   .. py:method:: add_n_particles(x, y, t, px, py, pt, qm, bunch_charge=None, w=None, sx=None, sy=None, sz=None)

      Add new particles to the container for fixed s.

      Either the charge (bunch_charge) of the particles added in this call, or the weight of each
      particle (w) must be provided.

      .. note::

         This can only be used *after* the initialization (grids) have
         been created, meaning after the call to :py:meth:`ImpactX.init_grids`
         has been made in the ImpactX class.

      .. attention::

         In MPI-parallel simulations, ``beam.add_n_particles(...)`` is local to the MPI rank, spatial locality does not matter.
         Thus, you can add particles at any MPI rank, e.g., equally chuncked up for perfect load balancing.

         You do NOT want to add the same unique particle at multiple MPI ranks.
         If you use the ``bunch_charge`` argument, then it will be interpreted as the charge of particles on the current rank.

         When ImpactX needs to sort particles spatially, it will redistribute them over MPI ranks automatically during tracking.

      :param x: positions in x
      :param y: positions in y
      :param t: positions as time-of-flight in c*t
      :param px: momentum in x
      :param py: momentum in y
      :param pt: momentum in t
      :param qm: charge over mass in 1/eV
      :param bunch_charge: total charge within a bunch in C
      :param w: weight of each particle: the macroparticle charge in units of the elementary charge `e` (i.e., how many real particles to represent)
      :param sx: spin component in x (optional; if provided, sy and sz must also be provided)
      :param sy: spin component in y (optional; if provided, sx and sz must also be provided)
      :param sz: spin component in z (optional; if provided, sx and sy must also be provided)

   .. py:method:: ref_particle()

      Access the reference particle (:py:class:`impactx.RefPart`).
      Deprecated: use :py:attr:`ref`.

      :return: return a data reference to the reference particle
      :rtype: impactx.RefPart

   .. py:property:: ref

      Access the reference particle (:py:class:`impactx.RefPart`).

   .. py:method:: set_ref_particle(refpart)

      Set reference particle attributes.

      :param impactx.RefPart refpart: a reference particle to copy all attributes from

   .. py:method:: set_bucket_length(bucket_length)

      Set length of the longitudinal particle domain (e.g., length of the RF bucket in z), optionally provided for the application of particle boundary conditions.

      :param bucket_length: length of the longitudinal particle domain in m

   .. py:method:: beam_moments()

      Calculate beam moments at current ``s`` like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters.

      :return: beam properties with string keywords
      :rtype: dict

   .. py:property:: store_beam_moments

      In situ calculate and store the beam moments for every simulation step (default: ``False``).

   .. py:method:: beam_moments_history()

      Return the history of the beam moments on every step.

      :return: Pandas Dataframe of beam properties, including the global reference position s
      :rtype: Pandas Dataframe

   .. py:method:: reset_beam_moments_history()

      Reset the history of the beam moments

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

   .. py:method:: to_df(local=True)

      Return all beam particles as a `pandas DataFrame <https://pandas.pydata.org/docs/user_guide/dsintro.html#dataframe>`__.

      Columns correspond to the particle attributes
      (``position_x``, ``position_y``, ``position_t``,
      ``momentum_x``, ``momentum_y``, ``momentum_t``,
      optional ``spin_x``, ``spin_y``, ``spin_z``,
      plus ``qm``, ``w``, and the AMReX-internal ``cpu`` and ``id``).
      Phase-space coordinates are relative to the reference particle (see :py:attr:`ref`).

      :param bool local: if ``True`` (default), only particles on the current MPI rank are returned.
                         If ``False``, particles are gathered to the root rank.
      :return: particle data as a ``pandas.DataFrame``, or ``None`` on ranks without particles.

      .. note::

         The returned DataFrame is a *copy* of the particle data; modifying it does not write back to the simulation.
         For in-place modification, iterate over particle tiles instead. See :ref:`usage-howto-python-particle-data` and
         the `pyAMReX particles guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__.

   .. py:attribute:: iterator

      Alias for :py:class:`impactx.ImpactXParIter`, the per-tile iterator used to access
      particle data on a mesh-refinement level without copying.
      The canonical usage is to import the iterator directly:

      .. code-block:: python

         from impactx import ImpactXParIter

         for pti in ImpactXParIter(sim.beam, level=0):
             soa = pti.soa().to_xp()  # NumPy (CPU) or CuPy (GPU)
             x = soa.real["position_x"]
             px = soa.real["momentum_x"]

      See :ref:`usage-howto-python-particle-data` for full usage examples and
      `pyAMReX particles <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__ for the underlying API.


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

   .. py:property:: qm_ratio_SI

      Read-only: Get reference particle charge to mass ratio (C/kg)

   .. py:method:: reset(keep_mass=False, keep_charge=False)

      Reset the reference particle.

      :param keep_mass: do not reset the reference particle mass
      :param keep_charge: do not reset the reference particle charge

   .. py:method:: set_species(species_name)

      Set reference particle species by name.
      Sets charge, mass, and gyromagnetic anomaly for a known particle species.
      Returns the reference particle for chaining.

      Known species: ``electron``, ``positron``, ``proton``, ``Hminus``.
      For other species, set charge, mass, and gyromagnetic anomaly individually via
      :py:meth:`set_charge_qe`, :py:meth:`set_mass_MeV`, and :py:meth:`set_gyromagnetic_anomaly`.

      .. dropdown:: Species Constants
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: ../../../src/particles/ReferenceParticle.H
            :language: cpp
            :dedent: 12
            :start-after: // [known_species]
            :end-before: // [/known_species]

      :param str species_name: particle species name

      Example usage:

      .. code-block:: python

         ref = sim.beam.ref
         ref.set_species("electron").set_kin_energy_MeV(2.0e3)

   .. py:method:: set_charge_qe(charge_qe)

      Write-only: Set reference particle charge in (positive) elementary charges.

   .. py:method:: set_mass_MeV(massE)

      Write-only: Set reference particle rest mass (MeV/c^2).

   .. py:method:: set_kin_energy_MeV(kin_energy_MeV)

      Write-only: Set reference particle kinetic energy (MeV)

   .. py:method:: load_file(madx_file)

      Load reference particle information from a MAD-X file.

      .. warning::

         Our MAD-X parser is under active development and provided as a preview.
         Please check any loaded MAD-X beams very carefully.
         Please report your experience and bugs on `our issue tracker <https://github.com/BLAST-ImpactX/impactx/issues>`__.

      :param madx_file: file name to MAD-X file with a ``BEAM`` entry


.. py:class:: impactx.ImpactXParIter(particle_container, level)

   Per-tile iterator over the beam particle data on a mesh-refinement level.
   Yields direct (zero-copy) access to the particle Struct-of-Arrays on CPU or GPU and is the
   recommended way to access particles in performance-critical paths.

   :param impactx.ParticleContainer particle_container: the beam particle container (e.g. ``sim.beam``).
   :param int level: mesh-refinement level (typically ``0``).

   See :ref:`usage-howto-python-particle-data` for full usage examples and the
   `pyAMReX particles guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__
   for the underlying API.

.. py:class:: impactx.ImpactXParConstIter(particle_container, level)

   Read-only variant of :py:class:`impactx.ImpactXParIter`.


Initial Beam Phase Space Distributions
--------------------------------------

This module provides particle beam distributions that can be used to initialize particle beams in an :py:class:`impactx.ParticleContainer`.

.. note::

    For additional information, consult the documentation on :ref:`theory-collective-beam-distribution-input`.
    For **all** except the ``thermal`` distribution we allow input in two forms:

    1. Phase space ellipse axis intersections (ImpactX native)
    2. Courant-Snyder (Twiss) parameters

For the input from Twiss parameters in Python, please use the helper function ``twiss``:

.. autofunction:: impactx.twiss

For computing Fourier coefficients from on-axis field data (used by :py:class:`~impactx.elements.RFCavity`, :py:class:`~impactx.elements.SoftQuadrupole`, and :py:class:`~impactx.elements.SoftSolenoid`):

.. autofunction:: impactx.fourier_coefficients

.. py:class:: impactx.distribution.Gaussian(lambdaX, lambdaY, lambdaT, lambdaPx, lambdaPy, lambdaPt, muxpx=0.0, muypy=0.0, mutpt=0.0, meanX=0.0, meanY=0.0, meanT=0.0, meanPx=0.0, meanPy=0.0, meanPt=0.0, dispX=0.0, dispPx=0.0, dispY=0.0, dispPy=0.0, cutX=0.0, cutY=0.0, cutT=0.0)

   A 6D Gaussian distribution, optionally with truncation.
   The user may specify an independent cutoff in each phase plane (x,px), (y,py), and (t,pt).
   The cut is performed in normalized Courant-Snyder variables corresponding to the user-supplied second moments or Twiss functions.
   As a result, this is equivalent to a cut corresponding to the (linearized) action in each plane.
   A cutoff value of 0 means no truncation (default).

   :param lambdaX: phase space position axis intercept; for zero correlation, these are the related RMS sizes (in meters)
   :param lambdaY: see lambdaX
   :param lambdaT: see lambdaX
   :param lambdaPx: phase space momentum axis intercept; for zero correlation, these are the related normalized RMS momenta (in radians)
   :param lambdaPy: see lambdaPx
   :param lambdaPt: see lambdaPx
   :param muxpx: correlation length-momentum
   :param muypy: see muxpx
   :param mutpt: see muxpx
   :param meanX: mean value of x-coordinate
   :param meanY: see meanX
   :param meanT: see meanX
   :param meanPx: mean value of x-momentum
   :param meanPy: see meanPx
   :param meanPt: see meanPx
   :param dispX: beam horizontal dispersion (in meters)
   :param dispPx: beam horizontal dispersion derivative (dimensionless)
   :param dispY: see dispX
   :param dispPy: see dispPx
   :param cutX: number of sigma at which to cut the distribution in (x,px) (dimensionless); 0 means no cut
   :param cutY: number of sigma at which to cut the distribution in (y,py) (dimensionless); 0 means no cut
   :param cutT: number of sigma at which to cut the distribution in (t,pt) (dimensionless); 0 means no cut

.. py:class:: impactx.distribution.Kurth4D(lambdaX, lambdaY, lambdaT, lambdaPx, lambdaPy, lambdaPt, muxpx=0.0, muypy=0.0, mutpt=0.0, meanX=0.0, meanY=0.0, meanT=0.0, meanPx=0.0, meanPy=0.0, meanPt=0.0, dispX=0.0, dispPx=0.0, dispY=0.0, dispPy=0.0)

   A 4D Kurth distribution transversely + a uniform distribution
   in t + a Gaussian distribution in pt.

.. py:class:: impactx.distribution.Kurth6D(lambdaX, lambdaY, lambdaT, lambdaPx, lambdaPy, lambdaPt, muxpx=0.0, muypy=0.0, mutpt=0.0, meanX=0.0, meanY=0.0, meanT=0.0, meanPx=0.0, meanPy=0.0, meanPt=0.0, dispX=0.0, dispPx=0.0, dispY=0.0, dispPy=0.0)

   A 6D Kurth distribution.

   R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)
   C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)

.. py:class:: impactx.distribution.KVdist(lambdaX, lambdaY, lambdaT, lambdaPx, lambdaPy, lambdaPt, muxpx=0.0, muypy=0.0, mutpt=0.0, meanX=0.0, meanY=0.0, meanT=0.0, meanPx=0.0, meanPy=0.0, meanPt=0.0, dispX=0.0, dispPx=0.0, dispY=0.0, dispPy=0.0)

   A K-V distribution transversely + a uniform distribution
   in t + a Gaussian distribution in pt.

.. py:class:: impactx.distribution.Empty

   This distribution sets all values to zero.

.. py:class:: impactx.distribution.Semigaussian(lambdaX, lambdaY, lambdaT, lambdaPx, lambdaPy, lambdaPt, muxpx=0.0, muypy=0.0, mutpt=0.0, meanX=0.0, meanY=0.0, meanT=0.0, meanPx=0.0, meanPy=0.0, meanPt=0.0, dispX=0.0, dispPx=0.0, dispY=0.0, dispPy=0.0)

   A 6D Semi-Gaussian distribution (uniform in position, Gaussian in momentum).

.. py:class:: impactx.distribution.Triangle(lambdaX, lambdaY, lambdaT, lambdaPx, lambdaPy, lambdaPt, muxpx=0.0, muypy=0.0, mutpt=0.0, meanX=0.0, meanY=0.0, meanT=0.0, meanPx=0.0, meanPy=0.0, meanPt=0.0, dispX=0.0, dispPx=0.0, dispY=0.0, dispPy=0.0)

   A triangle distribution for laser-plasma acceleration related applications.

   A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).
   The transverse distribution is a 4D waterbag.

.. py:class:: impactx.distribution.Waterbag(lambdaX, lambdaY, lambdaT, lambdaPx, lambdaPy, lambdaPt, muxpx=0.0, muypy=0.0, mutpt=0.0, meanX=0.0, meanY=0.0, meanT=0.0, meanPx=0.0, meanPy=0.0, meanPt=0.0, dispX=0.0, dispPx=0.0, dispY=0.0, dispPy=0.0)

   A 6D Waterbag distribution.

.. py:class:: impactx.distribution.Thermal(k, kT, kT_halo, normalize, normalize_halo, halo=0.0)

   A 6D stationary thermal or bithermal distribution.

Initial Beam Spin Distribution
------------------------------

.. py:class:: impactx.distribution.SpinvMF(mux, muy, muz)

   A von Mises-Fisher (vMF) distribution on the unit 2-sphere.

   This is used for initializing particle spin. There is a natural bijective correspondence between vMF distributions and mean (polarization) vectors.

   The algorithm used here is a simplification of the algorithm described in:
   C. Pinzon and K. Jung, "Fast Python sampler of the von Mises Fisher distribution", in the special case of the 2-sphere. Additional references used include:

   - K. V. Mardia and P. E. Jupp, Directional Statistics, Wiley, 1999;
   - S. Kang and H-S. Oh, "Novel sampling method for the von Mises-Fisher distribution", Stat. and Comput. 34, 106 (2024), `DOI:10.1007/s11222-024-10419-3 <https://doi.org/10.1007/s11222-024-10419-3>`__

   :param mux: x component of the unit vector specifying the mean direction
   :param muy: y component of the unit vector specifying the mean direction
   :param muz: z component of the unit vector specifying the mean direction


Lattice
-------

This module provides elements and methods for the accelerator lattice.

.. py:class:: impactx.elements.KnownElementsList

   An iterable, ``list``-like type of elements.

   .. py:method:: clear()

      Clear the list to become empty.

   .. py:method:: extend(list)

      Add a list of elements to the list.

   .. py:method:: append(element)

      Add a single element to the list.

   .. py:method:: load_file(filename, nslice=1, *, min_model="linear")

      Load and append a lattice file from MAD-X (.madx) or PALS (e.g., .pals.yaml) formats.

      The reader picks the *cheapest* ImpactX element that represents each imported
      element.
      ``min_model`` raises the lower gate of that choice, so that an imported lattice
      starts out with at least the requested level of physical fidelity
      (see :ref:`theory-assumptions`):

      * ``"linear"`` (default) uses the linear model
      * ``"paraxial"`` requires at least the chromatic ``Chr*`` models
      * ``"exact"`` requires at least the exact-Hamiltonian ``Exact*`` models

      A floor never *lowers* a model.
      An element that already requires a richer model, such as a thick octupole
      that only exists as ``ExactMultipole``, keeps it.
      Where a tier is not implemented for an element family, the next higher one is
      used.
      Where no model reaches the requested tier at all, the most faithful model
      available is used and a warning is emitted:

      .. list-table::
         :header-rows: 1
         :widths: 20 25 25 30

         * - MAD-X element
           - ``"linear"``
           - ``"paraxial"``
           - ``"exact"``
         * - ``DRIFT``
           - ``Drift``
           - ``ChrDrift``
           - ``ExactDrift``
         * - ``QUADRUPOLE``
           - ``Quad``
           - ``ChrQuad``
           - ``ExactQuad``
         * - ``SBEND``, ``RBEND``
           - ``Sbend`` / ``CFbend``
           - ``ExactSbend`` / ``ExactCFbend``
           - ``ExactSbend`` / ``ExactCFbend``
         * - ``DIPEDGE`` (and bend edges)
           - ``DipEdge(model="linear")``
           - ``DipEdge(model="nonlinear")``
           - ``DipEdge(model="nonlinear")``
         * - ``SOLENOID``
           - ``Sol``
           - ``ChrAcc(ez=0)``
           - ``ChrAcc(ez=0)`` (warns)
         * - ``SEXTUPOLE``, ``OCTUPOLE``, skew ``QUADRUPOLE``
           - ``ExactMultipole``
           - ``ExactMultipole``
           - ``ExactMultipole``

      A skew ``QUADRUPOLE`` reaches ``ExactMultipole`` at every floor because the
      reader translates a combined normal and skew quadrupole as one multipole,
      not because the cheaper tiers could not express it: a pure skew quadrupole
      is also a ``Quad`` or ``ChrQuad`` under a 45 degree ``rotation``.

      The paraxial ``SOLENOID`` is ``ChrAcc`` with ``ez=0``, which is the same
      hard-edge solenoid map for the reference particle and adds the chromatic
      dependence.
      ``Sol`` takes its strength per unit rigidity, as ``ks``, while ``ChrAcc``
      takes ``bz = ks * beta_gamma``, so the paraxial solenoid is pinned to the
      reference energy that the MAD-X ``BEAM`` command declares.
      Changing the reference energy after the import does not rescale it, and the
      reader warns when it applies this translation.
      Only an ``ENERGY`` that the file actually states is used for this: MAD-X
      defaults it to 1 GeV, so a lattice whose ``BEAM`` command omits it keeps
      ``Sol`` at every floor rather than converting against an assumed energy.

      Elements that MAD-X itself defines as thin (``MULTIPOLE``, ``KICKER``,
      ``RFCAVITY``, ``NLLENS``) and structural elements (``MARKER``, apertures,
      monitors) have no model tiers, so the kick or map they translate into is
      the same at every floor.
      Their *length*, however, is carried by drifts that the reader adds around
      that kick, and those follow ``min_model`` like any other drift.
      A finite-length ``RFCAVITY`` is a thin kick model, for example, becomes
      ``ExactDrift`` + ``ShortRF`` + ``ExactDrift`` at ``min_model="exact"``.
      The same applies to the drifts emitted for ``COLLIMATOR``, ``INSTRUMENT``,
      ``PLACEHOLDER`` and for a ``MONITOR`` with ``L > 0``.

      .. warning::

         Our MAD-X and PALS parsers are under active development and provided as a preview.
         Please check any loaded lattice files very carefully.
         Please report your experience and bugs on `our issue tracker <https://github.com/BLAST-ImpactX/impactx/issues>`__.

      .. seealso::

         *synmadx*, our :ref:`alternative MAD-X parser <usage-python-synmadx>` ported from Synergia.

         :py:meth:`impactx.elements.FilteredElementsList.replace_with_drifts` uses the
         same ``"linear"``/``"paraxial"``/``"exact"`` vocabulary for its ``model``
         parameter.

      :param filename: filename to file with beamline elements
      :param nslice: number of slices used for the application of collective effects
      :param min_model: lowest element model to translate into, one of ``"linear"``
         (default), ``"paraxial"`` or ``"exact"``

      **Example:**

      .. code-block:: python

         # import a lattice with at least exact-Hamiltonian element models
         sim.lattice.load_file("fodo.madx", nslice=25, min_model="exact")

   .. py:method:: from_pals(pals_line, nslice=1, *, min_model="linear")

      Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) Python Line.

      :param pals_line: PALS Python Line with beamline elements
      :param nslice: number of slices used for the application of collective effects
      :param min_model: lowest element model to translate into, see
         :py:meth:`~impactx.elements.KnownElementsList.load_file`

   .. py:method:: select(kind=None, name=None)

      Filter elements by type and/or name.
      If both are provided, OR-based logic is applied.

      Returns references to original elements, allowing modification and chaining.
      Chained ``.select(...).select(...)`` selections are AND-filtered.

      :param kind: Element type(s) to filter by. Can be a string (e.g., ``"Drift"``), regex pattern (e.g., ``r".*Quad"``), element type (e.g., ``elements.Drift``), or list/tuple of these.
      :param name: Element name(s) to filter by. Can be a string, regex pattern, or ``list``/``tuple`` of these.
      :rtype: :py:class:`impactx.elements.FilteredElementsList`

      **Examples:**

      .. code-block:: python

         # Filter by element type
         drift_elements = lattice.select(kind="Drift")
         quad_elements = lattice.select(kind=elements.Quad)

         # Filter by regex pattern
         all_quads = lattice.select(kind=r".*Quad")  # matches Quad, ChrQuad, ExactQuad

         # Filter by name
         specific_elements = lattice.select(name="quad1")

         # Chain filters (AND logic)
         drift_named_d1 = lattice.select(kind="Drift").select(name="drift1")

         # Modify original elements through references
         drift_elements[0].ds = 2.0  # modifies original lattice

         # delete all drifts
         lattice.select(kind=r".*Drift").delete()

         # replace all Quads with drift equivalents
         lattice.select(kind=r".*Quad").replace_with_drifts()

   .. py:method:: get_kinds()

      Get all unique element types in the lattice.

      :return: List of unique element types (sorted by name)
      :rtype: list[type]

   .. py:method:: count_by_kind(kind_pattern)

      Count elements of a specific kind.

      :param kind_pattern: Element kind to count. Can be string (e.g., "Drift"), regex pattern (e.g., r".*Quad"), or element type (e.g., elements.Drift)
      :return: Number of elements of the specified kind
      :rtype: int

   .. py:method:: has_kind(kind_pattern)

      Check if list contains elements of a specific kind.

      :param kind_pattern: Element kind to check for. Can be string (e.g., "Drift"), regex pattern (e.g., r".*Quad"), or element type (e.g., elements.Drift)
      :return: True if at least one element of the specified kind exists
      :rtype: bool

   .. py:method:: transfer_map(ref, order="linear", fallback_identity_map=False)

      Calculate the end-to-end transfer map of the elements in the list.

      Currently only the linear transfer map is implemented (``order="linear"``);
      the ``order`` parameter is reserved for future higher-order extensions.
      In linear mode the 6x6 map is composed element by element, using each
      element's analytic per-slice linear transport map.

      Collective effects like space charge, Coherent/Incoherent Synchrotron
      Radiation (CSR/ISR), and wakefield effects are not applied here; the
      returned map describes the purely linear single-particle dynamics of the
      design lattice.

      Phase-space ordering in the returned matrix is ``(x, px, y, py, t, pt)``.

      .. dropdown:: Example
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: tests/python/test_lattice_optics.py
            :language: bash

      :param ref: reference particle at the starting s
      :param order: So far, only the calculation of linear transfer maps is supported.
      :param fallback_identity_map: For elements with an undefined transfer map in the lattice, assume the identity matrix.
      :return: The end-to-end transfer map of the lattice.
      :rtype: Map6x6

   .. py:method:: map_trace(ref)

      Trace the cumulative 6x6 linear transport map element by element.

      The reference particle is passed by value (intentional copy); the
      caller's reference particle is not modified in place. This matches the
      convention used by :py:meth:`~transfer_map`.

      This per-element trace is what :py:meth:`~impactx.impactx_pybind.ImpactX.twiss`
      consumes to transport Twiss functions through the lattice.

      If you only need the final cumulative map at the lattice exit, prefer
      :py:meth:`~transfer_map` instead of indexing the last entry of
      :py:meth:`~map_trace`.

      :param ref: A reference particle.
      :return: A list of dictionaries, one per lattice element plus a leading
               entry for the starting position. Each entry contains:

               * ``s``    -- integrated path length along the reference
                 orbit, in meters;
               * ``name`` -- user-supplied element name (empty string if not
                 named);
               * ``type`` -- element type string (e.g. ``"Drift"``,
                 ``"Quad"``, ``"Sbend"``);
               * ``M``    -- cumulative 6x6 linear transport map from the
                 start of the lattice to the exit of this element (a
                 ``Map6x6`` instance; call ``.to_numpy()`` for a standard
                 C-ordered NumPy array).

               The first entry always has the identity map at the starting
               ``s``; the last entry contains the same map as
               :py:meth:`~transfer_map`.
      :rtype: list[dict]

   .. py:method:: to_dicts()

      Serialize the lattice to a list of dictionaries.

      Each element is converted to a dictionary using its ``to_dict()`` method.
      The resulting list can be serialized to JSON, YAML, or other formats.

      .. note::

         This transforms the buggy ``.to_dict()`` keys of
         ``ExactSbend``, ``PlaneXYRot``, ``PRot`` and ``ThinDipole``
         to degrees, which by accident are written in radians.
         See this comment in
         `issue #1367 <https://github.com/BLAST-ImpactX/impactx/issues/1367#issuecomment-4160236826>`__.

      :return: List of element dictionaries
      :rtype: list[dict]

      **Example:**

      .. code-block:: python

         import json
         from impactx import elements

         lattice = elements.KnownElementsList([
             elements.Drift(ds=1.0, name="d1"),
             elements.Quad(ds=0.5, k=2.0, name="q1"),
         ])

         # Serialize to JSON
         with open("lattice.impactx.json", "w") as f:
             json.dump(lattice.to_dicts(), f, indent=2)

   .. py:method:: from_dicts(dicts)

      Load and append elements from a list of dictionaries.

      Each dictionary should be in the format produced by ``to_dict()``,
      containing at minimum a ``type`` key identifying the element class.

      :param dicts: List of element dictionaries
      :type dicts: list[dict]

      **Example:**

      .. code-block:: python

         import json
         from impactx import elements

         # Load from JSON
         with open("lattice.impactx.json") as f:
             data = json.load(f)

         lattice = elements.KnownElementsList()
         lattice.from_dicts(data)

   .. py:method:: to_py()

      Generate Python code that recreates this lattice.

      Returns a string containing a complete Python script with imports
      and a ``get_lattice()`` function that returns a KnownElementsList
      with all elements.

      .. note::

         Like ``to_dicts()``, this transforms the buggy ``.to_dict()`` keys of
         ``ExactSbend``, ``PlaneXYRot``, ``PRot`` and ``ThinDipole``
         from radians to degrees.

      :return: Python source code
      :rtype: str

      **Example:**

      .. code-block:: python

         from impactx import elements

         lattice = elements.KnownElementsList([
             elements.Drift(ds=1.0, name="d1"),
             elements.Quad(ds=0.5, k=2.0, name="q1"),
         ])

         # Generate Python code
         code = lattice.to_py()
         print(code)

         # Save to file
         with open("my_lattice.py", "w") as f:
             f.write(code)

         # Later, use the generated file:
         # from my_lattice import get_lattice
         # lattice = get_lattice()

   .. py:method:: __eq__(other)

      Element-wise equality.
      Number of elements and every pair of elements must compare equal under ``==``.

   .. py:method:: isclose(other, *, rtol=1e-12, atol=0.0, ignore_attributes=None)

      Tolerant element-wise comparison. Number of elements must match.
      For each pair of elements at the same index, calls the element's own ``isclose``.

      :param other: Any iterable of elements
         (:py:class:`~impactx.elements.KnownElementsList`,
         :py:class:`~impactx.elements.FilteredElementsList`, or plain ``list``).
      :param rtol: Relative tolerance (default ``1e-12``).
      :param atol: Absolute tolerance (default ``0.0``).
      :param ignore_attributes: ``to_dict()`` keys to skip when comparing each
         pair of elements. Accepts a single string or any iterable of strings.
         Forwarded to each element's ``isclose``; see
         :ref:`element-comparison-methods` for the full semantics, including
         the special ``"type"`` key for cross-variant comparisons.
      :rtype: bool

   .. py:method:: plot_survey(ref=None, ax=None, legend=True, legend_ncols=5)

      Plot over s of all elements in the KnownElementsList.

      A positive element strength denotes horizontal focusing (e.g. for quadrupoles) and bending to the right (for dipoles).  In general, this depends on both the sign of the field and the sign of the charge.

      Either populates the matplotlib axes in ax or creates a new axes containing the plot.

      :param ref: A reference particle, checked for the charge sign to plot focusing/defocusing strength directions properly.
      :param ax: A plotting area in matplotlib (called axes there).
      :param legend: Plot a legend if true.
      :param legend_ncols: Number of columns for lattice element types in the legend.

.. py:class:: impactx.elements.FilteredElementsList

   View returned by :py:meth:`~impactx.elements.KnownElementsList.select` on a
   :py:class:`~impactx.elements.KnownElementsList` or by chained ``.select()`` on a filtered view.
   Indexing returns the same element objects as the full lattice; assigning to fields updates the
   underlying list.

   All mutating operations (``delete``, ``replace_each``, ``replace_with_drifts``) rebuild the
   lattice using cloned elements. Existing Python references to lattice elements will then point
   to objects that are **no longer in the lattice**. If you cache element references, re-fetch
   them from the lattice after any mutation.

   If the selection is empty, ``delete`` is a no-op and ``replace_*`` return
   an empty ``FilteredElementsList``.

   .. py:method:: select(kind=None, name=None)

      Narrow this view with an additional AND filter. OR logic within a single call matches
      :py:meth:`~impactx.elements.KnownElementsList.select`.

      :param kind: Same meaning as for :py:meth:`~impactx.elements.KnownElementsList.select`.
      :param name: Same meaning as for :py:meth:`~impactx.elements.KnownElementsList.select`.
      :rtype: :py:class:`impactx.elements.FilteredElementsList`

   .. py:method:: delete()

      Remove all elements in the current selection from the underlying lattice. Invalidates this
      view **and** all other live selections on that lattice.
      Call :py:meth:`~impactx.elements.KnownElementsList.select` on the underlying lattice again to obtain a new view.

      :rtype: None

   .. py:method:: replace_each(element, *, keep_name=True, keep_ds=False)

      Replace each selected element with a copy of ``element``. Invalidates all **other** live
      selections on the same lattice; returns a **new** filtered view over the same indices (the
      returned view is valid).

      :param element: Element to clone at each selected index (names and ``ds`` may be overridden;
         see below).
      :param keep_name: If true (default), copy ``name`` from each replaced element when present.
      :param keep_ds: If true, copy segment length ``ds`` from each replaced element; otherwise
         ``ds`` comes from the template (default false).
      :rtype: :py:class:`impactx.elements.FilteredElementsList`

      **Examples:**

      .. code-block:: python

         # Replace quadrupoles; names kept, ds and k from template
         sim.lattice.select(kind="Quad").replace_each(
             elements.Quad(name="tpl", ds=0.1, k=1.5)
         )

         # Same but keep ds from the replaced elements (only k from template)
         sim.lattice.select(kind="Quad").replace_each(
             elements.Quad(name="tpl", ds=0.1, k=1.5), keep_ds=True
         )

   .. py:method:: replace_with_drifts(*, model="match", keep_alignment=True, keep_aperture=False)

      Replace each selected element with a drift of the chosen physics family. Names and ``ds``
      are always taken from the replaced element. Invalidates all **other** live selections on the
      same lattice; returns a **new** filtered view over the same indices (the returned view is
      valid).

      :param model: With ``"match"`` (default), linear elements become ``Drift``, class names
         starting with ``Chr`` become ``ChrDrift``, and class names starting with ``Exact`` become
         ``ExactDrift``. With ``"linear"``, ``"paraxial"``, or ``"exact"``, every selected slot
         uses that drift type.
      :param keep_alignment: If true (default), copy ``dx``, ``dy``, and ``rotation`` from each
         replaced element; otherwise zero them.
      :param keep_aperture: If true, copy ``aperture_x`` and ``aperture_y`` from each replaced
         element; otherwise zero them (default false).
      :rtype: :py:class:`impactx.elements.FilteredElementsList`

      **Examples:**

      .. code-block:: python

         # All Quads become drifts of the matching model
         sim.lattice.select(kind=r".*Quad").replace_with_drifts()

         # Clear alignment errors and apertures
         sim.lattice.select(kind=r".*Quad").replace_with_drifts(keep_alignment=False)

   .. py:method:: __eq__(other)

      Element-wise equality, with the same semantics as
      :py:meth:`impactx.elements.KnownElementsList.__eq__`. A filtered view
      compares equal to any iterable of elements (plain Python ``list``,
      :py:class:`~impactx.elements.KnownElementsList`, or another
      :py:class:`~impactx.elements.FilteredElementsList`) as long as the
      lengths and pairwise element comparisons match.

   .. py:method:: isclose(other, *, rtol=1e-12, atol=0.0, ignore_attributes=None)

      Tolerant element-wise comparison, with the same semantics as
      :py:meth:`impactx.elements.KnownElementsList.isclose`.
      ``ignore_attributes`` is forwarded to each element's ``isclose``.

.. _element-comparison-methods:

Common comparison methods on lattice elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every lattice element class in :py:mod:`impactx.elements` supports the following value-based
comparison methods. They derive directly from each element's ``to_dict()`` output.

.. py:method:: impactx.elements.Element.__eq__(other)

   Value-based equality. Two elements are equal if they are instances of the
   same class and their ``to_dict()`` outputs match key-for-key. Float-valued
   fields are compared via ``==`` (so ``NaN != NaN``); list/matrix values are
   compared element-wise.

.. py:method:: impactx.elements.Element.isclose(other, *, rtol=1e-12, atol=0.0, ignore_attributes=None)

   Tolerant equality.
   Float-valued fields are compared via
   ``math.isclose(rel_tol=rtol, abs_tol=atol)``; lists/matrices of floats are compared
   element-wise. All other value types (ints, strings, ``None``) fall back to strict ``==``.
   Mismatched element types and foreign operands return ``False`` (unless ``"type"`` is
   listed in ``ignore_attributes``).

   :param other: Element to compare against.
   :param rtol: Relative tolerance forwarded to ``math.isclose`` / ``numpy.allclose``.
   :param atol: Absolute tolerance. Default ``0.0``.
   :param ignore_attributes: ``to_dict()`` keys to skip during the comparison.
      Accepts a single string or any iterable of strings. Useful when comparing
      loaded files where bookkeeping fields such as ``"name"`` should not
      affect the verdict. Including the special key ``"type"`` disables the
      same-class check, so e.g. :py:class:`~impactx.elements.Drift` and
      :py:class:`~impactx.elements.ExactDrift` can be compared on their common
      parameters; remaining keys must still match.
   :rtype: bool

   **Examples:**

   .. code-block:: python

      from impactx import elements

      a = elements.Drift(ds=1.0, name="d1")
      b = elements.Drift(ds=1.0 + 1e-15, name="d2")

      a == b      # False — name differs and ds differs by float noise
      a.isclose(b)                              # False — name still differs
      a.isclose(b, ignore_attributes="name")    # True

      # Compare a Drift to its exact-physics counterpart on common fields.
      lin = elements.Drift(ds=1.0)
      ex = elements.ExactDrift(ds=1.0)
      lin.isclose(ex, ignore_attributes=["type"])  # True

      # Lattice-wide comparison forwards ignore_attributes to each element.
      lattice_a.isclose(lattice_b, ignore_attributes=["name"])

.. py:function:: impactx.elements.isclose(a, b, *, rtol=1e-12, atol=0.0, ignore_attributes=None)

   Free-function form of :py:meth:`isclose`, equivalent to ``a.isclose(b, ...)``.
   Accepts either two elements or two iterables of elements
   (``KnownElementsList``, ``FilteredElementsList``, plain ``list``).

   .. code-block:: python

      from impactx import elements

      # two elements
      d1 = elements.Drift(ds=1.0, name="d1")
      d2 = elements.Drift(ds=1.0 + 1e-15, name="d2")
      elements.isclose(d1, d2, ignore_attributes="name")

      # two lattices (KnownElementsList, FilteredElementsList, or plain list)
      elements.isclose(lattice_a, lattice_b, ignore_attributes=["name"])


.. _usage-python-synmadx:

Synergia Migration Parser: synmadx
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ImpactX's primary MAD-X importer is
:py:meth:`impactx.elements.KnownElementsList.load_file`.  ImpactX also ships
*synmadx*, a custom compatibility parser ported from
`Synergia <https://github.com/fnalacceleratormodeling/synergia2>`__ for users
migrating existing Synergia workflows.

.. warning::

   The synmadx reader has `known correctness and coverage issues
   <https://github.com/BLAST-ImpactX/impactx/issues/1584>`__.  Use the primary
   :py:meth:`impactx.elements.KnownElementsList.load_file` importer for new and
   general-purpose MAD-X workflows, and validate converted lattices carefully
   when using synmadx for migration.

.. note::

   synmadx is only available if ImpactX was compiled with the CMake option
   ``-DImpactX_SYNMADX=ON``, which needs `Boost <https://www.boost.org>`__ installed.

.. py:class:: impactx.synmadx.MadX_reader

   Read and parse a MAD-X file into a Synergia lattice.

   .. py:method:: parse(string)

      Parse MAD-X input given as a string.

      :param string: MAD-X statements to parse

   .. py:method:: get_lattice(line_name)
                  get_lattice(line_name, filename)

      Return the named beamline as a Synergia lattice, optionally reading and
      parsing a MAD-X file first.

      :param line_name: name of the beamline (line or sequence) to extract
      :param filename: path to a MAD-X (``.madx``) file

.. py:function:: impactx.synmadx.syn2_to_impactx(lattice, init_monitor=True, final_monitor=True, order=Order.exact)

   Convert a Synergia lattice parsed with synmadx into a list of ImpactX
   lattice elements.

   :param lattice: Synergia lattice, e.g., from :py:meth:`impactx.synmadx.MadX_reader.get_lattice`; must include a reference particle
   :param init_monitor: prepend a ``BeamMonitor`` that records the initial beam
   :param final_monitor: append a ``BeamMonitor`` that records the final beam
   :param order: tracking model used for the converted elements: ``Order.linear``, ``Order.chr`` (chromatic), or ``Order.exact``

.. py:function:: impactx.synmadx.unroll_impactx_lattice(lattice)

   Return a Python source string that recreates the given list of ImpactX
   elements, e.g., for saving a converted lattice to a file.

   :param lattice: a list of ImpactX lattice elements

See the :ref:`Simple Booster (synmadx MAD-X Parser) example <examples-simple-booster-synmadx>`
for a complete workflow.


Lattice Elements
----------------

Lattice elements expose multiple methods, including: (i) ``push(pc)`` to advance particles ``pc`` (e.g., ``sim.beam``),
(ii) ``push(cm, ref)`` to advance a covariance matrix ``cm``, (iii) ``transfer_map(ref)`` that returns the
element's analytic 6x6 linear transport map (phase-space ordering ``(x, px, y, py, t, pt)``) for the reference
particle ``ref``, (iv) ``reverse()`` to reverse the element in place, and (v) ``to_dict()`` to serialize it.
For an element with ``nslice`` > 1, the pushes and maps refer to a single ``ds/nslice`` slice.

Many elements take a transverse aperture via ``aperture_x`` and ``aperture_y``: the ``Aperture``
collimator, and as a beam pipe most other elements.
They all follow the same convention. Each plane is bounded independently: a half-aperture of zero
or less removes the constraint in that plane only, while the other plane still cuts.
Bounding a single plane gives a jaw (slit) collimator, bounding both an iris; with only
``aperture_y`` set, a particle is lost when ``|y| > aperture_y`` at any ``x``, and the
``rectangular`` and ``elliptical`` shapes degenerate to the same slab.
The aperture is disabled entirely only if both planes are zero or less, which is the default.

.. py:class:: impactx.elements.CFbend(ds, rc, k, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

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
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.ConstF(ds, kx, ky, kt, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   A linear Constant Focusing element.

   :param ds: Segment length in m.
   :param kx: Focusing strength for x in 1/m.
   :param ky: Focusing strength for y in 1/m.
   :param kt: Focusing strength for t in 1/m.
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   .. py:property:: kx

      focusing x strength in 1/m

   .. py:property:: ky

      focusing y strength in 1/m

   .. py:property:: kt

      focusing t strength in 1/m

.. py:class:: impactx.elements.DipEdge(psi, rc, g, R=1, K0=pi**2/6, K1=0, K2=1, K3=1/6, K4=0, K5=0, K6=0, model="linear", location="entry", modify_ref_part=False, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   Edge focusing associated with bend entry or exit

   The model here is based on:

   * K. Hwang and S. Y. Lee, PRAB 18, 122401 (2015).

   as represented in the explicit, symplectic form provided in:

   * C. Mitchell and K. Hwang, in Proc. NAPAC2025, TUP040, Sacramento, CA (2025).

   Here, ``g`` denotes the magnetic gap, which is a length scale that sets the rate of decay of the fringe field.  The values ``K0`` - ``K6`` denote
   dimensionless field integrals, describing the shape of the fringe field, as defined in eqs. (28-34) of the first reference above.  In
   particular, ``K2`` is the well-known fringe field parameter denoted ``FINT`` in MAD-X.  The default values of the field integrals ``K0`` - ``K6`` are
   those given in eq. (52), corresponding to a ``tanh`` (i.e. logistic) field profile.

   When ``model = "linear"``, the linearized map is used.  This model is identical to:

   * K. L. Brown, SLAC Report No. 75 (1982)

   when expanded to first order in ``g/rc`` (gap / radius of curvature).

   By comparison, note that the MAD-X DIPEDGE element uses as input the half-gap ``HGAP = g/2``, and sets the default value ``FINT = 0`` (while
   the corresponding default value of ``K2`` is set to 1).

   Note that the nonlinear model includes a nonzero horizontal translation (depending on the field integral values) that is present even for a particle that begins on the ideal "hard-edge" reference
   trajectory.  For a beam, this will result in a centroid offset that will produce centroid oscillations in the downstream beamline. In practice, this can be avoided by aligning the downstream elements with
   the true horizontal position (after including the effect of the fringe field).  To model this correction, we allow two options in the dipedge model:

   * the option ``modify_ref_part = False`` (default), in which the shift due to the fringe field is applied to each beam particle phase space vector but not to the reference particle phase space vector --
   this model makes sense if the shift due to the fringe field is not considered in the baseline design, so that downstream elements are aligned with the "idealized" reference trajectory

   * the option ``modify_ref_part = True``, in which the shift due to the fringe field is applied to the reference particle phase space vector, but not to the beam particle phase space vector --
   this model makes sense if the shift due to the fringe field is considered as part of the baseline design, so that downstream elements are aligned with the "shifted" reference trajectory

   :param psi: Pole face angle [radians]
   :param rc: Radius of curvature [m]
   :param g: Gap parameter [m]
   :param R: Length scale used in fringe field integrals [m]
   :param K0: Fringe field integral [unitless]
   :param K1: Fringe field integral [unitless]
   :param K2: Fringe field integral [unitless]
   :param K3: Fringe field integral [unitless]
   :param K4: Fringe field integral [unitless]
   :param K5: Fringe field integral [unitless]
   :param K6: Fringe field integral [unitless]
   :param model: the fringe field model: ``linear`` (default) or ``nonlinear``
   :param location: the fringe field edge location: ``entry`` (default) or ``exit``
   :param modify_ref_part: apply fringe field to the reference particle ``True`` or ``False`` (default)
   :param dx: horizontal translation error [m]
   :param dy: vertical translation error [m]
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

.. py:class:: impactx.elements.Drift(ds, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   A drift.

   :param ds: Segment length in m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.ChrDrift(ds, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   A drift with chromatic effects included.  The Hamiltonian is expanded
   through second order in the transverse variables (x,px,y,py), with the exact pt
   dependence retained.

   :param ds: Segment length in m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.ExactDrift(ds, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   A drift using the exact nonlinear transfer map.

   :param ds: Segment length in m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.Kicker(xkick, ykick, unit="dimensionless", dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   A thin transverse kicker.

   :param xkick: horizontal kick strength (dimensionless OR T-m)
   :param ykick: vertical kick strength (dimensionless OR T-m)
   :param unit: specification of units (``"dimensionless"`` in units of the magnetic rigidity of the reference particle or ``"T-m"``)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

.. py:class:: impactx.elements.LinearMap(R, ds=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   A custom, linear transport matrix.

   The matrix elements :math:`R(i,j)` are indexed beginning with 1, so that :math:`i,j=1,2,3,4,5,6`.

   The matrix :math:`R` multiplies the phase space vector :math:`(x,px,y,py,t,pt)`, where coordinates :math:`(x,y,t)` have units of m
   and momenta :math:`(px,py,pt)` are dimensionless.  So, for example, :math:`R(1,1)` is dimensionless, and :math:`R(1,2)` has units of m.

   The internal tracking methods used by ImpactX are symplectic.  However, if a user-defined linear map :math:`R` is provided, it is
   up to the user to ensure that the matrix :math:`R` is symplectic.  Otherwise, this condition may be violated.

   :param R: a linear transport map to multiply with the phase space vector :math:`(x,px,y,py,t,pt)`.
   :param ds: length associated with a user-defined linear element (defaults to 0), in m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

   .. note::

      This element cannot be sliced, so a collective effect kick (space charge,
      CSR, ISR) is applied as ``K(ds/2) M(ds) K(ds/2)``: half a kick, the element
      map, then the other half, at two collective solves per element
      (see :py:attr:`~impactx.ImpactX.strang_split`).
      Model it as several shorter elements to resolve collective effects along
      its length.

.. py:class:: impactx.elements.Multipole(multipole, K_normal, K_skew, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   A general thin multipole element.

   :param multipole: index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
   :param K_normal: Integrated normal multipole coefficient (meter^(-m+1))
                    = ds * 1/(magnetic rigidity in T-m) * (derivative of order :math:`m-1` of :math:`B_y` with respect to :math:`x`)
   :param K_skew: Integrated skew multipole coefficient (meter^(-m+1))
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

.. py:class:: impactx.elements.ExactCFbend(ds, k_normal, k_skew, unit=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, int_order=2, mapsteps=10, nslice=1, name=None)

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

   :param ds: Segment length in m.
   :param k_normal: Array of normal multipole coefficients (in meter^(-m) OR in T/meter^(m-1) for m=1,2,3,..)
   :param k_skew: Array of skew multipole coefficients (in meter^(-m) OR in T/meter^(m-1) for m=1,2,3,...)
   :param unit: specification of units for multipole coefficients (by default, these are normalized by magnetic rigidity)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param int_order: the order used for symplectic integration (2, 4, or 6)
   :param mapsteps: number of integration steps per slice used for symplectic integration
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   The arrays k_normal and k_skew must have the same number of elements.

   .. py:property:: unit

      unit specification for multipole coefficients

   .. py:property:: int_order

      the order used for symplectic integration (2, 4, or 6)

   .. py:property:: mapsteps

      number of integration steps per slice used for symplectic integration


.. py:class:: impactx.elements.ExactMultipole(ds, k_normal, k_skew, unit=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, int_order=2, mapsteps=10, nslice=1, name=None)

   A thick Multipole magnet using the exact relativistic Hamiltonian, including all kinematic nonlinearities.
   The user must provide arrays containing normal and skew multipole coefficients, which can be specified up to arbitrarily high order.
   The fields are assumed to be uniform along the longitudinal beamline coordinate (hard-edge model).
   The coefficients must appear in the following sequence:

   dipole, quadrupole, sextupole, octupole, etc...

   (Note: Dipole coefficients are currently ignored, and will be supported in a separate combined-function dipole element.)

   Particle tracking is performed using symplectic integration based on the Hamiltonian splitting H = H_1 + H_2.
   Here H_1 is the nonlinear Hamiltonian for a drift (including the kinematic square root),
   and H_2 is the term containing the vector potential, which is a superposition of multipole contributions.

   :param ds: Segment length in m.
   :param k_normal: Array of normal multipole coefficients (in meter^(-m) OR in T/meter^(m-1) for m=1,2,3,..)
   :param k_skew: Array of skew multipole coefficients (in meter^(-m) OR in T/meter^(m-1) for m=1,2,3,...)
   :param unit: specification of units for multipole coefficients (by default, these are normalized by magnetic rigidity)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param int_order: the order used for symplectic integration (2, 4, or 6)
   :param mapsteps: number of integration steps per slice used for symplectic integration
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   .. py:property:: unit

      unit specification for multipole coefficients

   .. py:property:: int_order

      the order used for symplectic integration (2, 4, or 6)

   .. py:property:: mapsteps

      number of integration steps per slice used for symplectic integration

.. py:class:: impactx.elements.Empty

   This element does nothing.

.. py:class:: impactx.elements.NonlinearLens(knll, cnll, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

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
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

.. py:class:: impactx.elements.BeamMonitor(name, backend="default", encoding="g", period_sample_intervals=1, particles=True, beam_moments=True)

   A beam monitor, writing all beam particles at fixed ``s`` to openPMD files.

   If the same element ``name`` is used multiple times, then an output series is created with multiple outputs.

   The `I/O backend <https://openpmd-api.readthedocs.io/en/latest/backends/overview.html>`_ for `openPMD <https://www.openPMD.org>`_ data dumps.
   ``bp4``/``bp5`` is the `ADIOS2 I/O library <https://csmd.ornl.gov/adios>`_, ``h5`` is the `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_, and ``json`` is a `simple text format <https://en.wikipedia.org/wiki/JSON>`_.
   ``json`` only works with serial/single-rank jobs.
   By default, the first available backend in the order given above is taken.

   openPMD `iteration encoding <https://openpmd-api.readthedocs.io/en/latest/usage/concepts.html#iteration-and-series>`__ determines if multiple files are created for individual output steps or not.
   Variable based is an experimental encoding and needs the ADIOS2 backend to store more than one iteration.
   openPMD-api 0.17 and newer read it with random access (``Access.read_only``) as well as with linear access (``Access.read_linear``); readers older than 0.17 see only a single iteration.

   :param name: name of the series
   :param backend: I/O backend, e.g., ``bp``, ``h5``, ``json``
   :param encoding: openPMD iteration encoding: (v)ariable based, (f)ile based, (g)roup based (default)
   :param period_sample_intervals: for periodic lattice, only output every Nth period (turn)
   :param particles: output the individual particles of the beam (default: ``True``)
   :param beam_moments: output the beam moments (reduced beam characteristics) as openPMD attributes (default: ``True``)

   .. py:property:: name

      name of the series

   .. py:property:: particles

      Output the individual particles of the beam (default: ``True``).
      When set to ``False``, the beam monitor becomes a fast, small diagnostic that stores only per-output attributes (reference particle and, if enabled, beam moments).

   .. py:property:: beam_moments

      Output the beam moments (reduced beam characteristics) of the whole bunch as openPMD attributes on the ``beam`` particle species (default: ``True``).

   .. py:property:: nonlinear_lens_invariants

      Compute and output the invariants H and I within the nonlinear magnetic insert element.
      Requires ``particles=True``.

   .. py:property:: alpha

      Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
      Horizontal and vertical values must be equal.

   .. py:property:: beta

      Twiss beta (in meters) of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
      Horizontal and vertical values must be equal.

   .. py:property:: tn

      Dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I.

   .. py:property:: cn

      Scale factor (in meters^(1/2)) of the IOTA nonlinear magnetic insert element used for computing H and I.

.. py:class:: impactx.elements.Source(distribution, openpmd_path, active_once=True, load_ref_particle=True, load_step=None, load_step_index=None, name=None)

   A particle source.
   Currently, this only supports openPMD files from our :py:class:`impactx.elements.BeamMonitor`

   :param distribution: Distribution type of particles in the source. currently, only ``"openPMD"`` is supported
   :param openpmd_path: path to the openPMD series
   :param active_once: Inject particles only for the first lattice period. Default: ``True``
   :param load_ref_particle: Restore the reference particle from the species metadata of the openPMD file. Default: ``True``
   :param load_step: Which step (iteration) to load from the openPMD series, selected by step number. Default: ``None``
   :param load_step_index: Which step (iteration) to load from the openPMD series, selected by position in the file. Default: ``None``
   :param name: an optional name for the element

   .. note::

      ``load_step`` is the ImpactX step at which the :py:class:`impactx.elements.BeamMonitor`
      wrote the beam, which is stored as the openPMD iteration in the file.
      These step numbers are usually not consecutive, because the global step counter also advances in
      the elements between two monitors: list them with, e.g., `openpmd-ls <https://openpmd-api.readthedocs.io/en/0.17.1/utilities/cli.html>`__ before selecting one.
      A step that is not in the file is an error, which lists the steps that are in it.
      A negative value is an error, too: the step numbers in a file are not negative, use
      ``load_step_index`` to count back from the last step.

      ``load_step_index`` selects the step by position in the file instead: ``0`` is the first step
      and ``-1`` is the last, counting back from it as in Python, e.g. ``load_step_index=-2`` is the
      second to last step.
      Use this when the step numbers in the file are not known, e.g. ``load_step_index=-2`` to
      continue from the turn before the last one that a ring wrote.
      An index that reaches past either end of the file is an error, which lists the steps in it.

      Set at most one of ``load_step`` and ``load_step_index``.
      If neither is set, the last step in the file is loaded.

   .. note::

      With ``load_ref_particle``, the reference particle state (:math:`s`, positions, momenta, mass, charge
      and the gyromagnetic anomaly) is restored exactly from the iteration the particles are read from.
      Set ``load_ref_particle=False`` to load all particles relative to a manually configured reference
      particle: the relative particle coordinates in the file (in particular the momenta :math:`p_x`,
      :math:`p_y` and :math:`p_t`, which are normalized by the reference particle momentum) are then
      interpreted with respect to the existing reference particle.
      The existing reference particle in that case needs to be manually configured before the tracking loop.
      This option acts when particles are tracked or pushed (:py:func:`impactx.push`),
      it has no effect in envelope or reference-particle-only tracking.

   .. attention::

      In MPI-parallel simulations, reading the particles in the source file is distributed over all MPI ranks:
      each of the :math:`N` ranks reads a different contiguous chunk of :math:`1/N`-th of the particles, independent of their position.
      This balances the I/O and memory load between the ranks.

      When ImpactX needs to sort particles spatially, it will redistribute them over MPI ranks automatically during tracking.

.. py:class:: impactx.elements.Programmable(ds=0.0, nslice=1, name=None)

   A programmable beam optics element.

   This element can be programmed to receive callback hooks into Python functions.
   See :ref:`usage-howto-python-extend` for a worked example.

   :param ds: Segment length in m.
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   .. note::

      The Strang split (:py:attr:`~impactx.ImpactX.strang_split`) cannot halve this element's
      push, so it halves the collective effect kick (space charge, CSR, ISR)
      instead: ``K(ds/2) M(ds) K(ds/2)`` per slice, at two collective solves per
      slice.

   .. note::

      The ``Programmable`` element is intended for *replacing* a beamline element's
      particle push with custom Python code (e.g. a non-linear kick, tabulated map,
      or ML surrogate). The ``beam_particles`` callback is invoked per particle
      tile for performance, so the beam is presented in chunks, not as a single
      global container, which is a poor fit for observation/analysis.

      For **in-situ analysis** of the beam, prefer a :py:attr:`~impactx.ImpactX.hook`
      callback instead: ``sim.beam`` and methods like
      :py:meth:`~impactx.ParticleContainer.to_df` give you the full beam in one place.
      See :ref:`usage-howto-python-particle-data`.

   .. py:property:: push

      This is a function hook for pushing the whole particle container.
      Either this function is implemented or ``beam_particles`` and ``ref_particle`` are needed.
      This accepts a function or lambda with the following arguments:

      .. py:method:: user_defined_function(pc: impactx.ParticleContainer, step: int)

         This function is called for the particle container as it passes through the element.
         Note that the reference particle must be updated *before* the beam particles are pushed.

   .. py:property:: beam_particles

      This is a function hook for pushing all beam particles.
      This accepts a function or lambda with the following arguments:

      .. py:method:: user_defined_beam_function(pti: impactx.ImpactXParIter, refpart: impactx.RefPart)

         This function is called repeatedly for all particle tiles or boxes in the beam particle container.
         Particles can be pushed and are relative to the reference particle

   .. py:property:: ref_particle

      This is a function hook for pushing the reference particle.
      This accepts a function or lambda with the following argument:

      .. py:method:: user_defined_refpart_function(refpart: impactx.RefPart)

         This function is called for the reference particle as it passes through the element.
         The reference particle is updated *before* the beam particles are pushed.

.. py:class:: impactx.elements.Quad(ds, k, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   A Quadrupole magnet.

   :param ds: Segment length in m.
   :param k:  Quadrupole strength in m^(-2) (MADX convention)
              = (gradient in T/m) / (rigidity in T-m)
              k > 0 horizontal focusing
              k < 0 horizontal defocusing
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.ChrQuad(ds, k, unit=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   A Quadrupole magnet, with chromatic effects included.  The Hamiltonian is expanded
   through second order in the transverse variables (x,px,y,py), with the exact pt
   dependence retained.

   :param ds: Segment length in m.
   :param k:  Quadrupole strength in m^(-2) (MADX convention, if unit = 0)
              = (gradient in T/m) / (rigidity in T-m)
          OR  Quadrupole strength in T/m (MaryLie convention, if unit = 1)
              k > 0 horizontal focusing
              k < 0 horizontal defocusing
   :param unit: specification of units for quadrupole field strength
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   .. py:property:: k

      quadrupole strength in 1/m^2 (or T/m)

   .. py:property:: unit

      unit specification for quad strength

.. py:class:: impactx.elements.ExactQuad(ds, k, unit=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, int_order=2, mapsteps=10, nslice=1, name=None)

   A Quadrupole magnet using the exact relativistic Hamiltonian, including all kinematic nonlinearities.
   Particle tracking is performed using symplectic integration based on the Hamiltonian splitting H = H_1 + H_2.
   Here H_1 is the Hamiltonian for a linear quadrupole (containing all terms quadratic in the phase space variables),
   and H_2 is the remainder (including the kinematic square root).  This suggested splitting appears for example in:

   D. L. Bruhwiler et al, in Proc. of EPAC 98, pp. 1171-1173 (1998).
   E. Forest, J. Phys. A: Math. Gen. 39, 5321 (2006).

   :param ds: Segment length in m.
   :param k:  Quadrupole strength in m^(-2) (MADX convention, if unit = 0)
              = (gradient in T/m) / (rigidity in T-m)
          OR  Quadrupole strength in T/m (MaryLie convention, if unit = 1)
              k > 0 horizontal focusing
              k < 0 horizontal defocusing
   :param unit: specification of units for quadrupole field strength
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param int_order: the order used for symplectic integration (2, 4, or 6)
   :param mapsteps: number of integration steps per slice used for symplectic integration
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   .. py:property:: k

      quadrupole strength in 1/m^2 (or T/m)

   .. py:property:: unit

      unit specification for quad strength

   .. py:property:: int_order

      the order used for symplectic integration (2, 4, or 6)

   .. py:property:: mapsteps

      number of integration steps per slice used for symplectic integration

.. py:class:: impactx.elements.QuadEdge(k, unit=0, flag="entry", dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   Hard-edge nonlinear fringe field map for a Quadrupole.  This is a nonlinear symplectic map (derived from a third-order Lie generator), representing
   the effect of quadrupole entry or exit fringe fields in the hard-edge limit. This is an explicit symplectification of the Lie map that appears in eq
   (28) of:  E. Forest and J. Milutinovic, Nucl. Instrum. and Methods in Phys. Res. A 269, 474-482 (1988).

   :param k:  Quadrupole strength in m^(-2) (MADX convention, if unit = 0)
              = (gradient in T/m) / (rigidity in T-m)
          OR  Quadrupole strength in T/m (MaryLie convention, if unit = 1)
              k > 0 horizontal focusing
              k < 0 horizontal defocusing
   :param unit: specification of units for quadrupole field strength
   :param flag: location of edge (``"entry"`` or ``"exit"``)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

.. py:class:: impactx.elements.ChrPlasmaLens(ds, k, unit=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   An active cylindrically symmetric plasma lens, with chromatic effects included.
   The Hamiltonian is expanded through second order in the transverse variables
   (x,px,y,py), with the exact pt dependence retained.

   :param ds: Segment length in m.
   :param k:  focusing strength in m^(-2) (if unit = 0)
              = (azimuthal magnetic field gradient in T/m) / (rigidity in T-m)
              OR  azimuthal magnetic field gradient in T/m (if unit = 1)
   :param unit: specification of units for plasma lens focusing strength
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   .. py:property:: k

      plasma lens focusing strength in 1/m^2 (or T/m)

   .. py:property:: unit

      unit specification for plasma lens focusing strength

.. py:class:: impactx.elements.ChrAcc(ds, ez, bz, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

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
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

   .. py:property:: ez

      electric field strength in 1/m

   .. py:property:: bz

      magnetic field strength in 1/m

.. py:class:: impactx.elements.RFCavity(ds, escale, freq, phase, *, cos_coefficients=None, sin_coefficients=None, z=None, field_on_axis=None, ncoef=None, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, mapsteps=10, nslice=1, name=None)

   A radiofrequency cavity.  See :ref:`Models of Soft-Edge Elements <theory-softedge-elements>`.

   Provide **either** pre-computed Fourier coefficients (``cos_coefficients``, ``sin_coefficients``)
   **or** raw on-axis field data (``z``, ``field_on_axis``, ``ncoef``), not both.
   When the latter is given, Fourier coefficients are computed automatically
   using :func:`impactx.fourier_coefficients`.

   The units used for the on-axis longitudinal electric field are described in the documentation of ``escale`` below.  For example, if the values used to
   describe the on-axis electric field (as specified in ``cos_coefficients``, ``sin_coefficients``, or ``gradient_on_axis``) attain a peak on-axis value of 1, then the parameter
   ``escale``, which multiplies this profile, specifies the peak value of the longitudinal electric field gradient on-axis, divided by particle rest energy.

   In this case, ``escale`` has units of inverse meters.

   :param ds: Segment length in m.
   :param escale: scaling factor for on-axis RF electric field in 1/m
                  = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV)
   :param freq: RF frequency in Hz
   :param phase: RF driven phase in degrees
   :param cos_coefficients: array of ``float`` cosine coefficients in Fourier expansion of on-axis electric field Ez (optional); default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__
   :param sin_coefficients: array of ``float`` sine coefficients in Fourier expansion of on-axis electric field Ez (optional); default is a 9-cell TESLA superconducting cavity model from `DOI:10.1103/PhysRevSTAB.3.092001 <https://doi.org/10.1103/PhysRevSTAB.3.092001>`__
   :param z: array of longitudinal positions in m, covering the element from entry (``min(z)``) to exit (``max(z)``); the range is scaled to ``ds`` (alternative to Fourier coefficients)
   :param field_on_axis: array of on-axis electric field Ez values, typically normalized to a peak absolute value of 1; multiplied by ``escale`` (alternative to Fourier coefficients)
   :param ncoef: number of Fourier coefficients to compute (alternative to Fourier coefficients)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param mapsteps: number of integration steps per slice used for map and reference particle push in applied fields
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.Sbend(ds, rc, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   An ideal sector bend.

   :param ds: Segment length in m.
   :param rc: Radius of curvature in m.
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.ExactSbend(ds, phi, B=0.0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

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
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.Buncher(V, k, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   A short RF cavity element at zero crossing for bunching (MaryLie model).

   :param V: Normalized RF voltage drop V = Emax*L/(c*Brho)
   :param k: Wavenumber of RF in 1/m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

.. py:class:: impactx.elements.ShortRF(V, freq, phase=-90.0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   A short RF cavity element (MAD-X model).

   :param V: Normalized RF voltage V = maximum energy gain/(m*c^2)
   :param freq: RF frequency in Hz
   :param phase: RF synchronous phase in degrees (phase = 0 corresponds to maximum energy gain, phase = -90 corresponds go zero energy gain for bunching)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

.. py:class:: impactx.elements.SoftSolenoid(ds, bscale, *, cos_coefficients=None, sin_coefficients=None, z=None, field_on_axis=None, ncoef=None, unit=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, mapsteps=10, nslice=1, name=None)

   A soft-edge solenoid.  See :ref:`Models of Soft-Edge Elements <theory-softedge-elements>`.

   Provide **either** pre-computed Fourier coefficients (``cos_coefficients``, ``sin_coefficients``)
   **or** raw on-axis field data (``z``, ``field_on_axis``, ``ncoef``), not both.
   When the latter is given, Fourier coefficients are computed automatically
   using :func:`impactx.fourier_coefficients`.

   The units used for the on-axis longitudinal magnetic field data are determined by the parameter ``unit``.  For example, if the values used to
   describe the on-axis profile (as specified in ``cos_coefficients``, ``sin_coefficients``, or ``field_on_axis``) attain a peak on-axis value of 1, then the parameter
   ``bscale``, which multiplies this profile, specifies the peak value of the longitudinal magnetic field gradient on-axis.  If ``unit=0``, this is normalized by the magnetic rigidity.

   :param ds: Segment length in m.
   :param bscale: Scaling factor for on-axis magnetic field Bz in inverse meters (if unit = 0)
              = (magnetic field Bz in T) / (rigidity in T-m)
              OR  Magnetic field Bz in T (SI units, if unit = 1)
   :param cos_coefficients: array of ``float`` cosine coefficients in Fourier expansion of on-axis magnetic field Bz
            (optional); default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__
   :param sin_coefficients: array of ``float`` sine coefficients in Fourier expansion of on-axis magnetic field Bz
            (optional); default is a thin-shell model from `DOI:10.1016/J.NIMA.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__
   :param z: array of longitudinal positions in m, covering the element from entry (``min(z)``) to exit (``max(z)``); the range is scaled to ``ds`` (alternative to Fourier coefficients)
   :param field_on_axis: array of on-axis magnetic Bz field values, typically normalized to a peak absolute value of 1; multiplied by ``bscale`` (alternative to Fourier coefficients)
   :param ncoef: number of Fourier coefficients to compute (alternative to Fourier coefficients)
   :param unit: specification of units for scaling of the on-axis longitudinal magnetic field
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param mapsteps: number of integration steps per slice used for map and reference particle push in applied fields
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.Sol(ds, ks, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, nslice=1, name=None)

   An ideal hard-edge Solenoid magnet.

   :param ds: Segment length in m.
   :param ks: Solenoid strength in m^(-1) (MADX convention) in (magnetic field Bz in T) / (rigidity in T-m)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.PRot(phi_in, phi_out, name=None)

   Exact map for a pole-face rotation in the x-z plane.

   :param phi_in: angle of the reference particle with respect to the longitudinal (z) axis in the original frame in degrees
   :param phi_out: angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in degrees
   :param name: an optional name for the element

.. py:class:: impactx.elements.PlaneXYRot(angle, dx=0, dy=0, rotation=0, name=None)

   Map for a transverse rotation in the x-y plane (i.e., about the reference velocity vector).

   :param angle: nominal angle of rotation in the x-y plane, in degrees
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param name: an optional name for the element

.. py:class:: impactx.elements.Aperture(aperture_x, aperture_y, repeat_x=0, repeat_y=0, shift_odd_x=False, shape="rectangular", action="transmit", dx=0, dy=0, rotation=0, name=None)

   A thin collimator element, applying a transverse aperture boundary.

   A half-aperture of zero or less disables the boundary in that plane, the same
   convention as every other element's ``aperture_x``/``aperture_y`` above.
   Note that with ``action="absorb"`` a disabled plane makes the absorbing domain
   unbounded in that direction, so disabling both planes absorbs the whole beam,
   while for ``action="transmit"`` the element has no effect on the beam.

   :param aperture_x: horizontal half-aperture (rectangular or elliptical) in m; zero or less disables the horizontal boundary
   :param aperture_y: vertical half-aperture (rectangular or elliptical) in m; zero or less disables the vertical boundary
   :param repeat_x: horizontal period for repeated aperture masking (inactive by default) (meter)
   :param repeat_y: vertical period for repeated aperture masking (inactive by default) (meter)
   :param shift_odd_x: for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
   :param shape: aperture boundary shape: ``"rectangular"`` (default) or ``"elliptical"``
   :param action: aperture domain action: ``"transmit"`` (default) or ``"absorb"``
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param name: an optional name for the element

   .. py:property:: shape

      aperture type (rectangular, elliptical)

   .. py:property:: action

      aperture type (transmit, absorb)

   .. py:property:: aperture_x

      maximum horizontal coordinate in m; zero or less removes the horizontal boundary

   .. py:property:: aperture_y

      maximum vertical coordinate in m; zero or less removes the vertical boundary

.. py:class:: impactx.elements.PolygonAperture(vertices_x, vertices_y, min_radius2=0.0, repeat_x, repeat_y, shift_odd_x, action="transmit", dx=0, dy=0, rotation=0, name=None)

   This element defines a thin collimator element applying a transverse polygon aperture boundary defined by :math:`(x,y)` coordinates
   and optional radius below which all particles are transmitted. The vertices must define a closed curve and be ordered in the counter-clockwise direction.
   The first and last vertices must be identical.

   :param vertices_x: sequence of aperture boundary :math:`x` coordinates in m
   :param vertices_y: sequence of aperture boundary :math:`y` coordinates in m
   :param min_radius2: radius-squared of a circle fully inscribed by the polygon aperture (default 0) (meters-squared)
   :param repeat_x: horizontal period for repeated aperture masking (inactive by default) (meter)
   :param repeat_y: vertical period for repeated aperture masking (inactive by default) (meter)
   :param shift_odd_x: for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
   :param action: aperture domain action: ``"transmit"`` (default) or ``"absorb"``
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param name: an optional name for the element

   .. py:property:: min_radius2

      radius-squared of a fully inscribed circle. Particles with radius-squared less than this value are transmitted by the aperture and the polygon calculation is skipped.

      aperture type (transmit, absorb)

.. py:class:: impactx.elements.SoftQuadrupole(ds, gscale, *, cos_coefficients=None, sin_coefficients=None, z=None, gradient_on_axis=None, ncoef=None, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, mapsteps=10, nslice=1, name=None)

   A soft-edge quadrupole.  See :ref:`Models of Soft-Edge Elements <theory-softedge-elements>`.

   Provide **either** pre-computed Fourier coefficients (``cos_coefficients``, ``sin_coefficients``)
   **or** raw on-axis field/gradient data (``z``, ``gradient_on_axis``, ``ncoef``), not both.
   When the latter is given, Fourier coefficients are computed automatically
   using :func:`impactx.fourier_coefficients`.

   The units used for the on-axis quadrupole gradient are the same as those used for the quadrupole strength ``k`` in the element Quad.  For example, if the values used to
   describe the on-axis profile (as specified in ``cos_coefficients``, ``sin_coefficients``, or ``gradient_on_axis``) attain a peak on-axis value of 1, then the parameter
   ``gscale``, which multiplies this profile, specifies the peak value of the quadrupole field gradient on-axis, divided by the magnetic rigidity.

   In this case, ``gscale`` has units of inverse meters squared.

   :param ds: Segment length in m.
   :param gscale: Scaling factor for on-axis field gradient in inverse meters squared.
   :param cos_coefficients: array of ``float`` cosine coefficients in Fourier expansion of on-axis field gradient dBy/dx
            (optional); default is a tanh fringe field model based on `<http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__
   :param sin_coefficients: array of ``float`` sine coefficients in Fourier expansion of on-axis field gradient dBy/dx
            (optional); default is a tanh fringe field model based on `<http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf>`__
   :param z: array of longitudinal positions in m, covering the element from entry (``min(z)``) to exit (``max(z)``); the range is scaled to ``ds`` (alternative to Fourier coefficients)
   :param gradient_on_axis: array of on-axis field gradient dBy/dx values, typically normalized to a peak absolute value of 1; multiplied by ``gscale`` (alternative to Fourier coefficients)
   :param ncoef: number of Fourier coefficients to compute (alternative to Fourier coefficients)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param mapsteps: number of integration steps per slice used for map and reference particle push in applied fields
   :param nslice: number of slices used for the application of space charge
   :param name: an optional name for the element

.. py:class:: impactx.elements.SpinMap(v, A, ds=0, dx=0, dy=0, rotation=0, name=None)

   A custom, user-specified spin map that acts on the spin 3-vector :math:`(s_x,s_y,s_z)`.  Spin maps are specified in the Lie-algebraic form:

   .. math::

      \vec{s}_f = M(\zeta)\vec{s}_i,\quad\quad M(\zeta)=e^{v\cdot L}e^{A\Delta\zeta\cdot L}.

   Here :math:`v` is a 3-vector that defines the axis and angle of rotation at the phase space design point, and :math:`A` is a 3x6 matrix that defines the spin-orbit coupling for particles not on the design point.
   Also, :math:`\Delta\zeta=(x,p_x,y,p_y,t,p_t)` denotes the 6-vector of phase space variables as deviations from the design orbit. The quantities :math:`L_x`, :math:`L_y`, and :math:`L_z` are standard 3x3 matrices that define a basis for the Lie algebra :math:`so(3)`.

   The vector components :math:`v(i)` and the matrix elements :math:`A(i,j)` are indexed beginning with 1, so that :math:`i=1,2,3` and :math:`j=1,2,3,4,5,6`.
   The vector :math:`v` and the matrix :math:`A` are defaulted to zero, so only entries that differ from zero need to be specified.

   The matrix :math:`A` multiplies the phase space vector :math:`(x,p_x,y,p_y,t,p_t)`, where coordinates :math:`(x,y,t)` have units of m
   and momenta :math:`(p_x,p_y,p_t)` are dimensionless.  The three components output are dimensionless.  So, for example, :math:`A(1,1)` has units of 1/m, and :math:`A(1,2)` is dimensionless.
   All three components of :math:`v` are dimensionless.

   :param v: a 1-indexed, 3x1, axis-angle vector that defines the spin rotation at the phase space design point
   :param A: a 1-indexed, 3x6, spin-orbit coupling matrix to multiply with the phase space vector :math:`(x,p_x,y,p_y,t,p_t)` that defines the spin rotation for off-design particles
   :param ds: length associated with a user-defined linear element (defaults to 0), in m
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param name: an optional name for the element

   .. note::

      This element cannot be sliced, so a collective effect kick (space charge,
      CSR, ISR) is applied as ``K(ds/2) M(ds) K(ds/2)``: half a kick, the element
      map, then the other half, at two collective solves per element
      (see :py:attr:`~impactx.ImpactX.strang_split`).
      Model it as several shorter elements to resolve collective effects along
      its length.

.. py:class:: impactx.elements.ThinDipole(theta, rc, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   A general thin dipole element.

   :param theta: Bend angle (degrees)
   :param rc: Effective curvature radius (meters)
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

   Reference:

   * G. Ripken and F. Schmidt, Thin-Lens Formalism for Tracking, CERN/SL/95-12 (AP), 1995.

.. py:class:: impactx.elements.TaperedPL(k, taper, unit=0, dx=0, dy=0, rotation=0, aperture_x=0, aperture_y=0, name=None)

   A thin nonlinear plasma lens with transverse (horizontal) taper

   .. math::

      B_x = g \left( y + \frac{xy}{D_x} \right), \quad \quad B_y = -g \left(x + \frac{x^2 + y^2}{2 D_x} \right)

   where :math:`g` is the (linear) field gradient in T/m and :math:`D_x` is the targeted horizontal dispersion in m.

   :param k:  integrated focusing strength in m^(-1) (if unit = 0)
              = (length in m) * (magnetic field gradient :math:`g` in T/m) / (magnetic rigidity in T-m)
          OR  integrated focusing strength in T (if unit = 1)
              = (length in m) * (magnetic field gradient :math:`g` in T/m)
   :param taper: horizontal taper parameter in m^(-1)
              = 1 / (target horizontal dispersion :math:`D_x` in m)
   :param unit: specification of units for plasma lens focusing strength
   :param dx: horizontal translation error in m
   :param dy: vertical translation error in m
   :param rotation: rotation error in the transverse plane [degrees]
   :param aperture_x: horizontal half-aperture (elliptical) in m
   :param aperture_y: vertical half-aperture (elliptical) in m
   :param name: an optional name for the element

   .. py:property:: k

      integrated plasma lens focusing strength in 1/m (or T)

   .. py:property:: taper

      horizontal taper parameter in 1/m

   .. py:property:: unit

      unit specification for plasma lens focusing strength


Methods
"""""""

Each lattice element provides a ``.to_dict()`` method, which can be used to serialize its configuration.


.. py:function:: elements.transformation.insert_element_every_ds(list, ds, element)

   Insert an element every s into an element list.

   Splits up every element that is on s = N * ds for N>0.

   :param list: element lattice list
   :param ds: spacing in meters along s to add an element
   :param element: the extra element to add every s


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
