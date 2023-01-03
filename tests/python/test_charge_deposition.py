#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


import math

import matplotlib.pyplot as plt
import numpy as np

import amrex
import impactx


def test_charge_deposition(save_png=True):
    """
    Deposit charge and access/plot it
    """
    pp_amr = amrex.ParmParse("amr")
    pp_amr.addarr("n_cell", [16, 24, 32])

    sim = impactx.ImpactX()

    assert sim.n_cell == [16, 24, 32]

    sim.load_inputs_file("examples/fodo/input_fodo.in")
    sim.space_charge = True
    sim.slice_step_diagnostics = False

    # Future:
    # sim.ncell = [25, 25, 45]
    # sim.domain = amrex.RealBox([1., 2., 3.], [4., 5., 6.])
    print(f"sim.n_cell={sim.n_cell}")
    print(f"sim.domain={sim.domain}")

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()

    sim.evolve()

    rho = sim.rho(lev=0)
    rs = rho.sum_unique(comp=0, local=False)

    gm = sim.Geom(lev=0)
    dr = gm.data().CellSize()
    dV = np.prod(dr)

    beam_charge = dV * rs  # in C
    assert math.isclose(beam_charge, -1.0e-9, rel_tol=1.0e-8)

    half_z = sim.n_cell[2] // 2  # order: x,y,z

    f = plt.figure()
    ax = f.gca()
    ng = rho.nGrowVect
    for mfi in rho:
        bx = mfi.validbox()
        rbx = amrex.RealBox(bx, dr, gm.ProbLo())

        arr = rho.array(mfi)
        arr_np = np.array(arr, copy=False)  # indices: comp, z, y, x

        # shift box to zero-based local mfi index space
        half_z_local = half_z - bx.lo_vect[2]
        bx.shift(bx.lo_vect * -1)
        # check if the current tile contains the half-z plane
        if half_z_local < 0 or half_z_local > arr_np.shape[2]:
            continue

        comp = 0
        mu = 1.0e6  # m->mu
        im = ax.imshow(
            # arr_np[comp, half_z, ...] * dV,  # including guard
            arr_np[comp, half_z_local, ng[1] : -ng[1], ng[0] : -ng[0]]
            * dV,  # w/o guard
            origin="lower",
            aspect="auto",
            extent=[rbx.lo(0) * mu, rbx.hi(0) * mu, rbx.lo(1) * mu, rbx.hi(1) * mu],
        )
        cb = f.colorbar(im)
        cb.set_label(r"charge density  [C/m$^3$]")
        ax.set_xlabel(r"$x$  [$\mu$m]")
        ax.set_ylabel(r"$y$  [$\mu$m]")
        if save_png:
            plt.savefig("charge_deposition.png")
        else:
            plt.show()


# implement a direct script run mode, so we can run this directly too,
# with interactive matplotlib windows, w/o pytest
if __name__ == "__main__":
    test_charge_deposition(save_png=False)
