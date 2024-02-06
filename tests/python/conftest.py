# -*- coding: utf-8 -*-

import os

import pytest

import amrex.space3d as amr
import impactx

if impactx.Config.have_mpi:
    from mpi4py import MPI  # noqa

    print("loaded mpi4py")
else:
    print("NO mpi4py load")

# base path for input files
basepath = os.getcwd()


@pytest.fixture(autouse=True, scope="function")
def amrex_init(tmpdir):
    with tmpdir.as_cwd():
        amr.initialize(
            [
                # print AMReX status messages
                "amrex.verbose=2",
                # throw exceptions and create core dumps instead of
                # AMReX backtrace files: allows to attach to
                # debuggers
                "amrex.throw_exception=1",
                "amrex.signal_handling=0",
                # abort GPU runs if out-of-memory instead of swapping to host RAM
                "amrex.abort_on_out_of_gpu_memory=1",
                # do not rely on implicit host-device memory transfers
                "amrex.the_arena_is_managed=0",
            ]
        )
        yield
        if amr.initialized():
            amr.finalize()
