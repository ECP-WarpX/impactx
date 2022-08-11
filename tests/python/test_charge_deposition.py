# -*- coding: utf-8 -*-


import math

import amrex
import impactx
import matplotlib.pyplot as plt
import numpy as np


def test_charge_deposition():
    """
    Deposit charge and access/plot it
    """
    sim = impactx.ImpactX()

    sim.load_inputs_file("examples/fodo/input_fodo.in")
    sim.set_slice_step_diagnostics(False)

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()

    sim.evolve()

    rho = sim.rho(lev=0)
    rs = rho.sum_unique(comp=0, local=False)

    gm = sim.Geom(lev=0)
    dr = gm.data().CellSize()
    dV = np.prod(dr)

    beam_charge = dV*rs  # in C
    assert math.isclose(beam_charge, 1.0e-9)

    f = plt.figure()
    ax = f.gca()
    ng = rho.nGrowVect
    for mfi in rho:
        bx = mfi.validbox()
        rbx = amrex.RealBox(bx, dr, gm.ProbLo())

        arr = rho.array(mfi)
        arr_np = np.array(arr, copy=False)

        half_z = arr_np.shape[1] // 2  # indices: comp, z, y, x
        comp = 0
        mu = 1.e6  # m->mu
        im = ax.imshow(
            #arr_np[comp, half_z, ...] * dV,  # including guard
            arr_np[
                comp,
                half_z,
                ng[1]:-ng[1],
                ng[0]:-ng[0]] * dV,  # w/o guard
            origin='lower',
            aspect='auto',
            extent=[
                rbx.lo(0) * mu, rbx.hi(0) * mu,
                rbx.lo(1) * mu, rbx.hi(1) * mu
            ]
        )
        cb = f.colorbar(im)
        cb.set_label(r"charge density  [C/m$^3$]")
        ax.set_xlabel(r"$x$  [$\mu$m]")
        ax.set_ylabel(r"$y$  [$\mu$m]")
        plt.show()


# implement a direct script run mode, so we can run this directly too,
# with interactive matplotlib windows, w/o pytest
if __name__ == '__main__':
    test_charge_deposition()
