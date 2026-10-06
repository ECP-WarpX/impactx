#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import Config, ImpactX, distribution, elements

# This is the vertical-bending variant of the 2-bend dogleg lattice: every bend
# (and its dip-edge focusing) is rolled by 90 degrees, so the dispersion that
# the dogleg generates appears in the vertical (y) plane instead of horizontal
# (x). To obtain the x <-> y mirror image of the horizontal dogleg, the initial
# Twiss functions are exchanged x <-> y as well.
#
# The same lattice is run first with the envelope (covariance) solver and then
# with particle tracking, exercising both the linear-map and the particle path
# through rolled (tilted) elements.

kin_energy_MeV = 5.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles


def dogleg_vertical_distribution():
    """A 5 GeV electron Waterbag with x <-> y exchanged relative to the
    horizontal dogleg (see ``run_dogleg.py``)."""
    return distribution.Waterbag(
        lambdaX=1.3084093142e-5,  # exchanged with lambdaY
        lambdaY=2.2951017632e-5,  # exchanged with lambdaX
        lambdaT=5.5555553e-8,
        lambdaPx=2.803697378e-6,  # exchanged with lambdaPy
        lambdaPy=1.598353425e-6,  # exchanged with lambdaPx
        lambdaPt=2.000000000e-6,
        muxpx=0.933345606203060,  # == muypy, so the exchange leaves it unchanged
        muypy=0.933345606203060,
        mutpt=0.999999961419755,
    )


def dogleg_vertical_lattice():
    """The horizontal dogleg lattice with all bends and dip edges rolled by
    90 degrees, so it bends in the vertical (y) plane."""
    ns = 25  # number of slices per ds in the element
    rc = 10.3462283686195526  # bend radius (meters)
    psi = 0.048345620280243  # pole face rotation angle (radians)
    lb = 0.500194828041958  # bend arc length (meters)
    roll = 90.0  # roll the bending plane from horizontal (x) to vertical (y)

    # Drift elements
    dr1 = elements.Drift(name="dr1", ds=5.0058489435, nslice=ns)
    dr2 = elements.Drift(name="dr2", ds=0.5, nslice=ns)

    # Bend elements (rolled by 90 degrees)
    sbend1 = elements.Sbend(name="sbend1", ds=lb, rc=-rc, rotation=roll, nslice=ns)
    sbend2 = elements.Sbend(name="sbend2", ds=lb, rc=rc, rotation=roll, nslice=ns)

    # Dipole Edge Focusing elements (rolled by 90 degrees)
    dipedge1 = elements.DipEdge(
        name="dipedge1", psi=-psi, rc=-rc, g=0.0, K2=0.0, rotation=roll
    )
    dipedge2 = elements.DipEdge(
        name="dipedge2", psi=psi, rc=rc, g=0.0, K2=0.0, rotation=roll
    )

    return [sbend1, dipedge1, dr1, dipedge2, sbend2, dr2]


def assert_final_beam(moments):
    """Check the final beam against the dogleg reference values: the same
    final-beam properties that the horizontal dogleg checks in
    ``analysis_dogleg.py``, with x <-> y exchanged. The 90-degree roll makes the
    dogleg generate its dispersion in the vertical (y) plane, so it must show up
    in ``dispersion_y`` (and x <-> y throughout). Accepts the moments dict from
    either the envelope or the particle beam."""
    rtol = 2.0e-2
    if Config.precision == "SINGLE":
        # FIXME: sigma_y, sigma_t, emittance_y and emittance_t are off by 20% in
        # single precision, for both the envelope and the particle beam
        # https://github.com/BLAST-ImpactX/impactx/issues/1485
        rtol = 0.5
    assert np.allclose(
        [
            moments["sigma_x"],
            moments["sigma_y"],
            moments["sigma_t"],
            moments["emittance_x"],
            moments["emittance_y"],
            moments["emittance_t"],
            moments["alpha_x"],
            moments["beta_x"],
            moments["alpha_y"],
            moments["beta_y"],
            moments["dispersion_y"],
        ],
        [
            2.166654e-05,
            1.922660e-03,
            1.101353e-04,
            1.020439e-10,
            8.561046e-09,
            8.569865e-09,
            -1.347910e00,
            4.600362e00,
            1.338583e00,
            1.437407e01,
            -2.666972e-01,
        ],
        rtol=rtol,
        atol=0.0,  # do not let the default atol swallow the emittance checks
    )


# set up the simulation: reference particle, lattice, and diagnostics
sim = ImpactX()
sim.space_charge = False
sim.slice_step_diagnostics = True
sim.init_grids()

ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

monitor = elements.BeamMonitor("monitor", backend="h5")
sim.lattice.append(monitor)
sim.lattice.extend(dogleg_vertical_lattice())
sim.lattice.append(monitor)

# The same lattice is run with both solvers; init_envelope() copies the
# reference particle, so envelope tracking does not disturb the particle run.

# (1) envelope (covariance) tracking: exercises the linear-map (transport-map)
# path through the rolled elements.
sim.init_envelope(ref, dogleg_vertical_distribution(), bunch_charge_C)
sim.track_envelope()
assert_final_beam(sim.envelope.beam_moments(ref))

# (2) particle tracking: exercises the particle push through the rolled elements.
sim.add_particles(bunch_charge_C, dogleg_vertical_distribution(), npart)
sim.track_particles()
assert_final_beam(sim.beam.beam_moments())

sim.finalize()
