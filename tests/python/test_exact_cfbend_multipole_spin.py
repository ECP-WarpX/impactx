#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


import numpy as np

from impactx import Config, ImpactX, distribution, elements


def test_exact_cfbend_multipole_spin():
    sim = ImpactX()

    # set numerical parameters and IO control
    sim.particle_shape = 2  # B-spline order
    sim.space_charge = False
    sim.spin = True
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    # basic beam parameters
    kin_energy_MeV = 2.0e3
    bunch_charge_C = 1.0e-9  # used with space charge
    npart = 10000  # number of macro particles

    # set reference particle
    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
    distr = distribution.Waterbag(
        lambdaX=4.0e-5,
        lambdaY=4.0e-5,
        lambdaT=1.0e-3,
        lambdaPx=3.0e-5,
        lambdaPy=3.0e-5,
        lambdaPt=2.0e-4,
        muxpx=0.0,
        muypy=0.0,
        mutpt=0.0,
    )
    spin = distribution.SpinvMF(
        0.4,
        0.9,
        0.1,
    )

    sim.add_particles(bunch_charge_C, distr, npart, spin)

    # initial beam
    beam = sim.beam
    initial_beam_df = beam.to_df()

    # design the accelerator lattice
    ns = 1  # number of slices per ds in the element

    # add beam diagnostics
    monitor = elements.BeamMonitor("monitor", backend="h5")

    # lattice elements
    cfbend1 = elements.ExactCFbend(
        name="cfbend1",
        ds=1.0,
        k_normal=[5.0e-8, 1.0, -2.0],
        k_skew=[0.0, -0.5, 1.4],
        unit=0,
        int_order=2,
        mapsteps=20,
        nslice=ns,
    )

    multipole1 = elements.ExactMultipole(
        name="quad1",
        ds=-1.0,
        k_normal=[0.0, 1.0, -2.0],
        k_skew=[0.0, -0.5, 1.4],
        unit=0,
        int_order=2,
        mapsteps=20,
        nslice=ns,
    )

    h = 0.1
    L = 1.0
    cfbend2 = elements.ExactCFbend(
        name="cfbend2",
        ds=L,
        k_normal=[h],
        k_skew=[0.0],
        unit=0,
        int_order=2,
        mapsteps=20,
        nslice=ns,
    )

    phi_val = 180.0 * h * L / np.pi
    sbend1 = elements.ExactSbend(
        name="sbend1",
        ds=-L,
        phi=-phi_val,
        nslice=ns,
    )

    # set the lattice
    sim.lattice.append(monitor)
    sim.lattice.append(cfbend1)
    sim.lattice.append(multipole1)
    sim.lattice.append(cfbend2)
    sim.lattice.append(sbend1)

    # run simulation
    sim.track_particles()

    # initial beam
    beam = sim.beam
    final_beam_df = beam.to_df()

    # clean shutdown
    sim.finalize()

    # analysis
    PHASE_COLS = [
        "position_x",
        "position_y",
        "position_t",
        "momentum_x",
        "momentum_y",
        "momentum_t",
    ]
    SPIN_COLS = ["spin_x", "spin_y", "spin_z"]

    phase_atol = 1.0e-7
    if Config.precision == "SINGLE":
        # FIXME: the roundtrip does not close in single precision, position_t is off
        # by up to 1.6e-3
        # https://github.com/BLAST-ImpactX/impactx/issues/1483
        phase_atol = 5.0e-3
    for c in PHASE_COLS:
        np.testing.assert_allclose(
            final_beam_df[c].to_numpy(),
            initial_beam_df[c].to_numpy(),
            atol=phase_atol,
            rtol=0,
            err_msg=f"Roundtrip mismatch in {c}",
        )
    spin_atol = 3.0e-7
    if Config.precision == "SINGLE":
        # FIXME: the roundtrip does not close in single precision, spin_z is off by
        # up to 2.1e-4
        # https://github.com/BLAST-ImpactX/impactx/issues/1483
        spin_atol = 1.0e-3
    for c in SPIN_COLS:
        np.testing.assert_allclose(
            final_beam_df[c].to_numpy(),
            initial_beam_df[c].to_numpy(),
            atol=spin_atol,
            rtol=0,
            err_msg=f"Roundtrip mismatch in {c}",
        )
