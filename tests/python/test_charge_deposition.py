# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX
import matplotlib.pyplot as plt
import numpy as np


def test_charge_deposition():
    """
    Deposit charge and access/plot it
    """
    impactX = ImpactX()

    impactX.load_inputs_file("examples/fodo/input_fodo.in")

    impactX.init_grids()
    impactX.init_beam_distribution_from_inputs()
    impactX.init_lattice_elements_from_inputs()

    impactX.evolve()

    f = plt.figure()
    ax = f.gca()
    lev = 0
    rho = impactX.rho(lev)
    print(f"rho={rho}")
    print(f"rho.sum={rho.sum(0)}")

    gm = impactX.Geom(lev)
    print(f"gm={gm}")

    for mfi in rho:
        bx = mfi.validbox()
        #rbx = amrex.RealBox(bx, gm.CellSize(), gm.ProbLo())
        print(f"bx={bx}")
        #print(f"rbx={rbx}")

        arr = rho.array(mfi)
        arr_np = np.array(arr, copy=False)

        half_z = arr_np.shape[1] // 2
        im = ax.imshow(
            arr_np[0, half_z, ...]
            #extent=rbx
        )
        f.colorbar(im)
        plt.show()


if __name__ == '__main__':
    amrex.initialize([
        # print AMReX status messages
        "amrex.verbose=2",
        # throw exceptions and create core dumps instead of
        # AMReX backtrace files: allows to attach to
        # debuggers
        "amrex.throw_exception=1",
        "amrex.signal_handling=0",
        # abort GPU runs if out-of-memory instead of swapping to host RAM
        "amrex.abort_on_out_of_gpu_memory=1"
    ])

    test_charge_deposition()

    amrex.finalize()
