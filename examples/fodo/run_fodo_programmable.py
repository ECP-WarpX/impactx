#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

try:
    import cupy as cp

    cupy_available = True
except ImportError:
    cupy_available = False

import amrex
from impactx import Config, ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=3.9984884770e-5,
    sigmaY=3.9984884770e-5,
    sigmaT=1.0e-3,
    sigmaPx=2.6623538760e-5,
    sigmaPy=2.6623538760e-5,
    sigmaPt=2.0e-3,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# number of slices per ds in each lattice element
ns = 25


# build a custom, Pythonic beam optical element
def my_drift(pge, pti, refpart):
    """This pushes the beam particles as a drift.

    Relative to the reference particle.

    :param pti: particle iterator for the current tile or box
    :param refpart: the reference particle
    """
    # CPU/GPU logic
    if Config.have_gpu:
        if cupy_available:
            array = cp.array
        else:
            print("Warning: GPU found but cupy not available! Try managed...")
            array = np.array
        if Config.gpu_backend == "SYCL":
            print("Warning: SYCL GPU backend not yet implemented for Python")

    else:
        array = np.array
    # access AoS data such as positions and cpu/id
    aos = pti.aos()
    aos_arr = array(aos, copy=False)

    # access SoA data such as momentum
    soa = pti.soa()
    real_arrays = soa.GetRealData()
    px = array(real_arrays[0], copy=False)
    py = array(real_arrays[1], copy=False)
    pt = array(real_arrays[2], copy=False)

    # length of the current slice
    slice_ds = pge.ds / pge.nslice

    # access reference particle values to find beta*gamma^2
    pt_ref = refpart.pt
    betgam2 = pt_ref**2 - 1.0

    # advance position and momentum (drift)
    aos_arr[:]["x"] += slice_ds * px[:]
    aos_arr[:]["y"] += slice_ds * py[:]
    aos_arr[:]["z"] += (slice_ds / betgam2) * pt[:]


def my_ref_drift(pge, refpart):
    """This pushes the reference particle.

    :param refpart: reference particle
    """
    #  assign input reference particle values
    x = refpart.x
    px = refpart.px
    y = refpart.y
    py = refpart.py
    z = refpart.z
    pz = refpart.pz
    t = refpart.t
    pt = refpart.pt
    s = refpart.s

    # length of the current slice
    slice_ds = pge.ds / pge.nslice

    # assign intermediate parameter
    step = slice_ds / (pt**2 - 1.0) ** 0.5

    # advance position and momentum (drift)
    refpart.x = x + step * px
    refpart.y = y + step * py
    refpart.z = z + step * pz
    refpart.t = t - step * pt

    # advance integrated path length
    refpart.s = s + slice_ds


pge1 = elements.Programmable()
pge1.nslice = ns
pge1.beam_particles = lambda pti, refpart: my_drift(pge1, pti, refpart)
pge1.ref_particle = lambda refpart: my_ref_drift(pge1, refpart)
pge1.ds = 0.25
pge1.threadsafe = True  # allow OpenMP threading for speed

# attention: assignment is a reference for pge2 = pge1

pge2 = elements.Programmable()
pge2.nslice = ns
pge2.beam_particles = lambda pti, refpart: my_drift(pge2, pti, refpart)
pge2.ref_particle = lambda refpart: my_ref_drift(pge2, refpart)
pge2.ds = 0.5
pge2.threadsafe = True  # allow OpenMP threading for speed

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", "h5")

# design the accelerator lattice
fodo = [
    monitor,
    pge1,  # equivalent to elements.Drift(ds=0.25, nslice=ns)
    monitor,
    elements.Quad(ds=1.0, k=1.0, nslice=ns),
    monitor,
    pge2,  # equivalent to elements.Drift(ds=0.5, nslice=ns)
    monitor,
    elements.Quad(ds=1.0, k=-1.0, nslice=ns),
    monitor,
    pge1,  # equivalent to elements.Drift(ds=0.25, nslice=ns)
    monitor,
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
