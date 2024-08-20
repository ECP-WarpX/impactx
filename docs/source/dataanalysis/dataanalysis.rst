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
* ``t_ref`` clock time * c in meters
* ``px_ref`` momentum in x, normalized to mass*c, :math:`p_x = \gamma \beta_x`
* ``py_ref`` momentum in y, normalized to mass*c, :math:`p_y = \gamma \beta_y`
* ``pz_ref`` momentum in z, normalized to mass*c, :math:`p_z = \gamma \beta_z`
* ``pt_ref`` energy, normalized by rest energy, :math:`p_t = -\gamma`
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


.. _dataanalysis-beam-characteristics:

Reduced Beam Characteristics
----------------------------

ImpactX calculates reduced beam characteristics based on the beam moments during runtime.
These include averaged positions, momenta, beam emittances and Courant-Snyder (Twiss) parameters.
For computing beam moments (as elsewhere), positions and momenta are given as deviations with respect to the reference particle (see :ref:`Coordinates and Units <theory-coordinates-and-units>`).

The reduced beam characteristics are stored with the output of the beam monitor element.
They are also calculated before, after, and during each step of the simulation.
If ``diag.slice_step_diagnostics`` is enabled, they will also be calculated during each slice of each beamline element.

The code writes out the values in an ASCII file prefixed ``reduced_beam_characteristics`` containing the follow columns:

* ``step``
    Iteration within the simulation
* ``s``
    Reference particle coordinate ``s`` (unit: meter)
* ``x_mean/min/max``, ``y_mean/min/max``, ``t_mean/min/max``
    Average / minimum / maximum particle displacement with respect to the reference particle in the dimensions of ``x``, ``y`` (transverse coordinates, unit: meter), and ``t`` (normalized time difference :math:`ct`, unit: meter)
* ``sig_x``, ``sig_y``, ``sig_t``
    Standard deviation of the particle positions (speed of light times time delay for ``t``) (unit: meter)
* ``px_mean/min/max``, ``py_mean/min/max``, ``pt_mean/min/max``
    Average / minimum / maximum particle momentum deviation from the reference particle momentum, divided by the magnitude of the reference particle momentum (unit: dimensionless, radians for transverse momenta)
* ``sig_px``, ``sig_py``, ``sig_pt``
    Standard deviation of the particle momentum deviations (energy difference for ``pt``) normalized by the magnitude of the reference particle momentum (unit: dimensionless)
* ``emittance_x``, ``emittance_y``, ``emittance_t``
    Normalized rms beam emittance (unit: meter)
* ``alpha_x``, ``alpha_y``, ``alpha_t``
    Courant-Snyder (Twiss) alpha (unit: dimensionless).  Transverse Twiss functions are calculated after removing correlations with particle energy.
* ``beta_x``, ``beta_y``, ``beta_t``
    Courant-Snyder (Twiss) beta (unit: meter).  Transverse Twiss functions are calculated after removing correlations with particle energy.
* ``dispersion_x``, ``dispersion_y``
    Horizontal and vertical dispersion (unit: meter)
* ``dispersion_px``, ``dispersion_py``
    Derivative of horizontal and vertical dispersion (unit: dimensionless)
* ``charge``
    Total beam charge (unit: Coulomb)


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
