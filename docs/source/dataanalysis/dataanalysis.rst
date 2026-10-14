.. _dataanalysis:

Data Analysis
=============

.. _dataanalysis-monitor:

Beam Monitor
------------

ImpactX provides a zero-sized beam monitor element that can be placed in lattices to output the particle beam at multiple positions in a lattice.
Output is written in the standardized, `open particle-mesh data schema (openPMD) <https://www.openPMD.org>`__ and is `compatible with many codes and data analysis frameworks <https://github.com/openPMD/openPMD-projects>`__.

For data analysis of openPMD data, see examples of `many supported tools, Python libraries and frameworks <https://openpmd-api.readthedocs.io/en/latest/analysis/viewer.html>`__.
`Exporting data to ASCII <https://openpmd-api.readthedocs.io/en/latest/analysis/pandas.html#openpmd-to-ascii>`__ is possible, too.

See also `WarpX' documentation on openPMD <https://warpx.readthedocs.io/en/latest/dataanalysis/formats.html>`__.

At each monitor, the output for every particle in the beam is provided.
This includes the 6 canonical phase space variables (x [m], px, y[m], py, t[m], pt), where each coordinate is measured relative to the reference particle, in the local moving frame.
See the section :ref:`Coordinates and Units <theory-coordinates-and-units>` for additional details.

.. _dataanalysis-monitor-fixed-t:

Particles at Fixed t
"""""""""""""""""""""

With :pp:param:`<monitor_name>.coord_system = t`, a monitor writes the particles at fixed :math:`t` instead: a snapshot of the beam at the time at which the reference particle passes the monitor.
Each particle is drifted to that time.
Codes that use time as the independent variable can read this output, e.g., `WarpX <https://warpx.readthedocs.io>`__, and so can a ``source`` element with :pp:param:`<source_name>.coord_system = t`.

The 6 phase space variables are then (x [m], px, y [m], py, z [m], pz), in the units of the openPMD standard:

* ``position``: ``x``, ``y``, ``z`` in meters, relative to the ``positionOffset``, which is the position ``x_ref``, ``y_ref``, ``z_ref`` of the reference particle.
* ``momentum``: ``x``, ``y``, ``z`` in units of the reference momentum :math:`p_0 = \beta\gamma m c`.
  Their openPMD attribute ``unitSI`` is :math:`p_0` in kg·m/s.
  ``z`` is the total longitudinal momentum :math:`p_z/p_0` and not a deviation from the reference momentum.
* The openPMD ``time`` of the output is the time of the reference particle, in seconds.

The axes are those of the local frame that moves with the reference particle: ``z`` points along the reference momentum.

All other particle attributes, the reference particle and the bunch properties below are the same as at fixed :math:`s`.
In particular, the bunch properties describe the beam at fixed :math:`s`.

.. _dataanalysis-monitor-refparticle:

Additional Beam Attributes
""""""""""""""""""""""""""

We add the following additional attributes on the openPMD ``beam`` species at the monitor position.

Reference particle:

* ``beta_ref`` reference particle normalized velocity :math:`\beta = v/c`
* ``gamma_ref`` reference particle Lorentz factor :math:`\gamma = 1/\sqrt{1-\beta^2}`
* ``beta_gamma_ref`` reference particle momentum normalized to rest mass :math:`\beta\gamma = p/(mc)`
* ``s_ref`` integrated orbit path length, in meters
* ``x_ref`` horizontal position x, in meters
* ``y_ref`` vertical position y, in meters
* ``z_ref`` longitudinal position z, in meters
* ``t_ref`` clock time * c, in meters
* ``px_ref`` momentum in x, normalized to mass*c, :math:`p_x = \gamma \beta_x`
* ``py_ref`` momentum in y, normalized to mass*c, :math:`p_y = \gamma \beta_y`
* ``pz_ref`` momentum in z, normalized to mass*c, :math:`p_z = \gamma \beta_z`
* ``pt_ref`` energy, normalized by rest energy, :math:`p_t = -\gamma`
* ``gyromagnetic_anomaly_ref`` `anomalous magnetic moment <https://en.wikipedia.org/wiki/Anomalous_magnetic_dipole_moment>`__, unitless
* ``mass_ref`` reference rest mass, in kg
* ``charge_ref`` reference charge, in C

Bunch properties: all properties listed in :ref:`Reduced Beam Characteristics <dataanalysis-beam-characteristics>`.

Example to print the integrated orbit path length ``s`` at each beam monitor position:

.. code-block:: python

   import openpmd_api as io

   series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

   for k_i, i in series.iterations.items():
       beam = i.particles["beam"]
       s_ref = beam.get_attribute("s_ref")
       print(f"step {k_i:>3}: s_ref={s_ref}")


.. _dataanalysis-python-inmemory:

In-Memory Python Access
-----------------------

In addition to the file-based openPMD output above, ImpactX can be steered from a Python script
and the beam particles, the reference particle, and the space-charge fields can be inspected or
modified **live** during a run. This is the right approach for custom in-situ diagnostics,
turn-by-turn feedback, AI/ML coupling, or any analysis that you don't want to round-trip through (slow) filesystem disks.

* :ref:`usage-howto-python-extend`: the overall pattern, plus the callback surfaces
  (``sim.hook[...]`` and the :py:class:`~impactx.elements.Programmable` element).
  For analysis, prefer ``sim.hook[...]`` with ``sim.beam``: the ``Programmable``
  element hands you particles in per-tile chunks and is intended for custom
  particle pushes, not observation.
* :ref:`usage-howto-python-particle-data`: accessing ``sim.beam`` as a pandas DataFrame
  (via :py:meth:`~impactx.ParticleContainer.to_df`) or iterating over particle tiles for
  high-performance access.
* :ref:`usage-howto-python-field-data`: accessing the space-charge fields
  (:py:meth:`~impactx.ImpactX.rho`, :py:meth:`~impactx.ImpactX.phi`,
  :py:meth:`~impactx.ImpactX.space_charge_field`) and their lifetime caveats.


.. _dataanalysis-beam-characteristics:

Reduced Beam Characteristics
----------------------------

ImpactX calculates reduced beam characteristics based on the beam moments during runtime.
These include averaged positions and momenta, as well as beam emittances and Courant-Snyder (Twiss) parameters.
For computing beam moments (as elsewhere), positions and momenta are given as deviations with respect to the reference particle (see :ref:`Coordinates and Units <theory-coordinates-and-units>`).

The reduced beam characteristics are stored with the output of the beam monitor element.
They are also calculated before, after, and during each step of the simulation.
If :pp:param:`diag.slice_step_diagnostics` is enabled, they will also be calculated during each slice of each beamline element.

The code writes out the values in an ASCII file prefixed ``reduced_beam_characteristics`` containing the follow columns:

* ``step``
    Iteration within the simulation
* ``period``
    Current period (turn or cycle) for periodic lattices (unit: dimensionless)
* ``s``
    Reference particle coordinate ``s`` (unit: meter)
* ``mean/min/max_x``, ``mean/min/max_y``, ``mean/min/max_t``
    Average / minimum / maximum particle displacement with respect to the reference particle in the dimensions of ``x``, ``y`` (transverse coordinates, unit: meter), and ``t`` (normalized time difference :math:`ct`, unit: meter)
* ``sigma_x``, ``sigma_y``, ``sigma_t``
    Standard deviation of the particle positions (speed of light times time delay for ``t``) (unit: meter)
* ``mean/min/max_px``, ``mean/min/max_py``, ``mean/min/max_pt``
    Average / minimum / maximum particle momentum deviation from the reference particle momentum, divided by the magnitude of the reference particle momentum (unit: dimensionless, radians for transverse momenta)
* ``sigma_px``, ``sigma_py``, ``sigma_pt``
    Standard deviation of the particle momentum deviations (energy difference for ``pt``) normalized by the magnitude of the reference particle momentum (unit: dimensionless)
* ``emittance_x``, ``emittance_y``, ``emittance_t``
    Unnormalized rms beam emittances (unit: meter)
* ``alpha_x``, ``alpha_y``, ``alpha_t``
    Courant-Snyder (Twiss) alpha (unit: dimensionless).  Transverse Twiss functions are calculated after removing correlations with particle energy.
* ``beta_x``, ``beta_y``, ``beta_t``
    Courant-Snyder (Twiss) beta (unit: meter).  Transverse Twiss functions are calculated after removing correlations with particle energy.
* ``dispersion_x``, ``dispersion_y``
    Horizontal and vertical dispersion (unit: meter)
* ``dispersion_px``, ``dispersion_py``
    Derivative of horizontal and vertical dispersion (unit: dimensionless)
* ``emittance_xn``, ``emittance_yn``, ``emittance_tn``
    Normalized rms beam emittances (unit: meter)
* ``emittance_1``, ``emittance_2``, ``emittance_3``
    Normalized rms beam eigenemittances (aka mode emittances) (unit: meter)
    These three diagnostics are written optionally if the flag eigenemittances = True.
* ``charge_C``
    Total beam charge (unit: Coulomb)
* ``mean_sx``, ``mean_sy``, ``mean_sz``
    Mean of the particle spin vector (unit: dimensionless).  This is also known as the polarization vector.  This diagnostic is written only if spin tracking is on (algo.spin = True).
* ``sigma_sx``, ``sigma_sy``, ``sigma_sz``
    Standard deviation of the three spin components (unit: dimensionless).  This diagnostic is written only if spin tracking is on (algo.spin = True).

Quantities that a beam does not define are written as ``nan``.
The ``step``, ``period`` and ``s`` columns describe the position in the lattice and are always written.
The Courant-Snyder (Twiss) ``alpha`` and ``beta`` of a plane are ratios to its rms emittance, and are undefined wherever that emittance vanishes: for a single particle, for a cold plane, or longitudinally for a beam without energy spread.
For zero-weight test beams, added with a bunch charge of zero, there is no total weight to normalize by: every weighted moment is undefined and ``charge_C`` is zero, while the unweighted ``min`` and ``max`` still report the extent of the beam.
For a simulation whose beam is loaded in a :py:class:`~impactx.elements.Source` element instead of prior to tracking, the beam is empty until that element is reached: every diagnostic before it reports ``nan`` for all moments, minima and maxima included, and zero for ``charge_C``.


.. _dataanalysis-plot:

Interactive Analysis
--------------------

When steering ImpactX from Python, one can at any point visualize the beam phase space with:

.. code-block:: python

   import matplotlib.pyplot as plt

   from impactx import ImpactX, RefPart, distribution, elements

   sim = ImpactX()

   # ... setup and simulate ...

   pc = sim.particle_container()

   fig = pc.plot_phasespace()

   # note: figure data available on MPI rank zero
   if fig is not None:
       fig.savefig("phase_space.png")
       plt.show()

.. figure:: https://user-images.githubusercontent.com/1353258/295041638-8410ba76-9bd2-4dae-9810-5ec9f33dd372.png
   :alt: In situ visualization of the beam phase space projections.

   In situ visualization of the beam phase space projections.
