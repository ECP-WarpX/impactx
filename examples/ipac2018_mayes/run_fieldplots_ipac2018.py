#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Christopher E. Mayes
# License: BSD-3-Clause-LBNL
#
# Benchmark as in Figure 1 of
#   C. E. Mayes, R. D. Ryne, and D. C. Sagan, "3D Space Charge in BMAD," in proc. IPAC2018,
#   DOI:10.18429/JACoW-IPAC2018-THPAK085
#
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import amrex.space3d as amr
from impactx import ImpactX, RefPart, distribution, elements

# theory data from eq. (1) in Mayes, sampled to text files in
txt = np.loadtxt("Ex_Mayes.dat")
theory = pd.read_csv(
    "Ex_Mayes.dat", sep=" ", names=["x", "r=0.01", "r=0.1", "r=1", "r=10"]
)
x_theory = theory["x"]


def run_and_plot(case, r, ax):
    """Run ImpactX for a new r and plot

    Parameters
    ----------
    case : string
      case name
    r : Float
      aspect ratio = sigma_z / sigma_perp:  range 0.01 to 10
    ax : matplotlib.axes
      a figure axis to plot into
    """
    sim = ImpactX()

    # set numerical parameters and IO control
    sim.n_cell = [128, 128, 256]  # 128^3 or higher for reasonable convergence
    # sim.n_cell = [32, 32, 64]
    sim.particle_shape = 2  # B-spline order
    sim.space_charge = True
    sim.dynamic_size = True
    sim.prob_relative = 3.0

    # space charge solver fine-tuning
    sim.mlmg_relative_tolerance = 1.0e-8
    sim.mlmg_max_iters = 500
    sim.mlmg_verbosity = 2

    # beam diagnostics
    sim.diagnostics = False
    sim.slice_step_diagnostics = False

    # domain decomposition & space charge mesh
    sim.init_grids()

    # load a cold 10 MeV electron beam
    mass_MeV = 0.510998950  # electron mass in MeV/c^2
    energy_MeV = 1.00001 * mass_MeV  # reference energy (total)
    bunch_charge_C = 1.0e-9  # charge in C
    npart = int(1000000)  # number of macro particles

    #   reference particle (Mayes Fig. 1 is for + charge)
    ref = sim.particle_container().ref_particle()
    kinetic_MeV = energy_MeV - mass_MeV
    ref.set_charge_qe(1.0).set_mass_MeV(mass_MeV).set_energy_MeV(kinetic_MeV)

    #   particle bunch
    sigma_r = 1.0e-3  # fixed at 1 mm
    sigma_z = r * sigma_r
    gamma = energy_MeV / mass_MeV
    beta = (1.0 - (1.0 / gamma) ** 2) ** 0.5
    print(gamma, beta)
    c0 = 2.99792458e8  # speed of light in m/s
    sigma_t = sigma_z / (beta * gamma)  # recall t is implicitly scaled by c0
    print(f"sigma_t={sigma_t}m")

    distr = distribution.Gaussian(
        sigmaX=sigma_r,
        sigmaY=sigma_r,
        sigmaT=sigma_t,
        sigmaPx=0.0,
        sigmaPy=0.0,
        sigmaPt=0.0,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    # add beam diagnostics
    monitor = elements.BeamMonitor("monitor", backend="h5")

    # design the accelerator lattice
    sim.lattice.extend([monitor])

    # run simulation
    sim.evolve()

    # plotting

    # simulation data
    F_x = sim.space_charge_field(lev=0, comp="x")  # [V/m] in the lab frame

    gm = sim.Geom(lev=0)
    dr = gm.data().CellSize()
    dV = np.prod(dr)

    half_x, half_y, half_z = [n // 2 for n in sim.n_cell]  # order: x,y,z

    # plot data slices
    ng = F_x.nGrowVect
    q_e = -1.602176634e-19
    for mfi in F_x:
        bx = mfi.validbox()
        rbx = amr.RealBox(bx, dr, gm.ProbLo())

        arr_np = F_x.array(mfi).to_numpy(copy=True)  # indices: x, y, z, comp

        # crop off guard cells
        arr_np = arr_np[ng[0] : -ng[0], ng[1] : -ng[1], ng[2] : -ng[2]]

        # shift box to zero-based local mfi index space
        half_y_local = half_y - bx.lo_vect[1]
        half_z_local = half_z - bx.lo_vect[2]

        # check if the current tile contains the half-y plane
        if half_y_local < 0 or half_y >= bx.hi_vect[1]:
            continue
        # check if the current tile contains the half-z plane
        if half_z_local < 0 or half_z >= bx.hi_vect[2]:
            continue

        comp = 0
        lineout = arr_np[..., half_y_local, half_z_local, comp]  # w/o guard, N
        E_x = lineout / gamma / 1.0e6  # E_x(x[m]) [MV/m]
        ax.plot(
            np.linspace(rbx.lo(0) / sigma_r, rbx.hi(0) / sigma_r, E_x.shape[0]),
            E_x,
            "ko",
            markersize=3,
        )
    xmin = -6.0
    xmax = 6.0
    ymin = -5.0
    ymax = 5.0
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.plot(x_theory, theory[case], label=case)  # E_x(sigma_r) [MV/m]

    # clean shutdown
    del sim


f = plt.figure()
ax = f.gca()

# run_and_plot("r=0.01", 0.01, ax)
run_and_plot("r=0.1", 0.1, ax)
run_and_plot("r=1", 1, ax)
run_and_plot("r=10", 10, ax)

# cb = f.colorbar(im)
# cb.set_label(r"charge density  [C/m$^3$]")
ax.set_xlabel(r"$x$  [$\sigma_x$]")
ax.set_ylabel(r"$E_x$  [MV/m]")
ax.legend()
# plt.savefig("charge_deposition.png")
plt.show()

# clean shutdown
amr.finalize()
