.. _dataanalysis:

Data Analysis
=============

Beam Monitor
------------

ImpactX provides a zero-sized beam monitor element that can be placed in lattices to output the particle beam at multiple positions in a lattice.
Output is written in the standardized, `open particle-mesh data schema (openPMD) <https://www.openPMD.org>`__ and is `compatible with many codes and data analysis frameworks <https://github.com/openPMD/openPMD-projects>`__.

For data analysis of openPMD data, see examples of `many supported tools, Python libraries and frameworks <https://openpmd-api.readthedocs.io/en/latest/analysis/viewer.html>`__.
`Exporting data to ASCII <https://openpmd-api.readthedocs.io/en/latest/analysis/pandas.html#openpmd-to-ascii>`__ is possible, too.

See also `WarpX' documentation on openPMD <https://warpx.readthedocs.io/en/latest/dataanalysis/formats.html>`__.

Additional Beam Attributes
""""""""""""""""""""""""""

We add the following additional attributes on the openPMD ``beam`` species at the monitor position.

Reference particle:

* ``beta_ref`` reference particle normalized velocity :math:`\beta = v/c`
* ``gamma_ref`` reference particle Lorentz factor :math:`\gamma = 1/\sqrt{1-\beta^2}`
* ``s_ref`` integrated orbit path length, in meters
* ``x_ref`` horizontal position x, in meters
* ``y_ref`` vertical position y, in meters
* ``z_ref`` longitudinal position z, in meters
* ``t_ref`` clock time * c in meters
* ``px_ref`` momentum in x, normalized to proper velocity :math:`p_x = \gamma v_x`
* ``py_ref`` momentum in y, normalized to proper velocity :math:`p_y = \gamma v_y`
* ``pz_ref`` momentum in z, normalized to proper velocity :math:`p_z = \gamma v_z`
* ``pt_ref`` energy deviation, normalized by rest energy
* ``mass`` reference rest mass, in kg
* ``charge`` reference charge, in C

Example to print the integrated orbit path length ``s`` at each beam monitor position:

.. code-block:: python

   import openpmd_api as io

   series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

   for k_i, i in series.iterations.items():
       beam = i.particles["beam"]
       s_ref = beam.get_attribute("s_ref")
       print(f"step {k_i:>3}: s_ref={s_ref}")


Reduced Beam Characteristics
----------------------------

ImpactX calculates reduced beam characteristics like averaged positions, momenta, beam emittances and Courant-Snyder (Twiss) parameters during runtime.
These quantities are calculated before, after, and during each step of the simulation.
If ``diag.slice_step_diagnostics`` is enabled, they will also be calculated during each slice of each beamline element.

The code writes out the values in an ASCII file prefixed ``reduced_beam_characteristics`` containing the follow columns:

* ``step``
    Iteration within the simulation
* ``s``, ``ref_beta_gamma``
    Reference particle coordinate ``s`` and relativistic momentum normalized by the particle mass and the speed of light
* ``x_mean``, ``y_mean``, ``t_mean``
    Average beam particle position in the dimensions of ``x``, ``y`` (transverse coordinates), and ``t`` (normalized time difference)
* ``sig_x``, ``sig_y``, ``sig_t``
    RMS of the average beam particle positions
* ``px_mean``, ``py_mean``, ``pt_mean``
    Average beam momenta
* ``sig_px``, ``sig_py``, ``sig_pt``
    RMS of the average beam momenta (energy difference for ``pt``)
* ``emittance_x``, ``emittance_y``, ``emittance_t``
    Normalized beam emittance
* ``alpha_x``, ``alpha_y``, ``alpha_t``
    Twiss alpha
* ``beta_x``, ``beta_y``, ``beta_t``
    Twiss beta
* ``charge``
    Cumulated beam charge in C
