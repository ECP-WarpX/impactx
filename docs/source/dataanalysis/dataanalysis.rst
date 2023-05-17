.. _dataanalysis:

Data Analysis
=============

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
    Cumulated beam charge

.. note::

   openPMD beam monitor output: TODO :-)

   For further information, please see https://warpx.readthedocs.io/en/latest/dataanalysis/formats.html.
