# -*- coding: utf-8 -*-

import amrex
import impactx
import pytest

if impactx.Config.have_mpi:
    from mpi4py import MPI
    print("loaded mpi4py")
else:
    print("NO mpi4py load")

@pytest.fixture(autouse=True, scope='session')
def amrex_init():
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
    yield
    amrex.finalize()
