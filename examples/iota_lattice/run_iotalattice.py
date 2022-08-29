#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.set_particle_shape(2)  # B-spline order
sim.set_space_charge(False)
# sim.set_diagnostics(False)  # benchmarking
sim.set_slice_step_diagnostics(True)

# domain decomposition & space charge mesh
sim.init_grids()

# init particle beam
energy_MeV = 2.5
charge_C = 1.0e-9  # assign zero weighting to particles
mass_MeV = 938.27208816
qm_qeeV = 1.0e-6 / mass_MeV
npart = 10000

distr = distribution.Waterbag(
    sigmaX=1.588960728035e-3,
    sigmaY=2.496625268437e-3,
    sigmaT=1.0e-3,
    sigmaPx=2.8320397837724e-3,
    sigmaPy=1.802433091137e-3,
    sigmaPt=0.0,
)
sim.add_particles(qm_qeeV, charge_C, distr, npart)

# set the energy in the reference particle
sim.particle_container().ref_particle().set_energy_MeV(energy_MeV, mass_MeV)

# init accelerator lattice
ns = 10  # number of slices per ds in the element

# Drift elements
dra1 = elements.Drift(ds=0.9125, nslice=ns)
dra2 = elements.Drift(ds=0.135, nslice=ns)
dra3 = elements.Drift(ds=0.725, nslice=ns)
dra4 = elements.Drift(ds=0.145, nslice=ns)
dra5 = elements.Drift(ds=0.3405, nslice=ns)
drb1 = elements.Drift(ds=0.3205, nslice=ns)
drb2 = elements.Drift(ds=0.14, nslice=ns)
drb3 = elements.Drift(ds=0.1525, nslice=ns)
drb4 = elements.Drift(ds=0.31437095, nslice=ns)
drc1 = elements.Drift(ds=0.42437095, nslice=ns)
drc2 = elements.Drift(ds=0.355, nslice=ns)
dnll = elements.Drift(ds=1.8, nslice=ns)
drd1 = elements.Drift(ds=0.62437095, nslice=ns)
drd2 = elements.Drift(ds=0.42, nslice=ns)
drd3 = elements.Drift(ds=1.625, nslice=ns)
drd4 = elements.Drift(ds=0.6305, nslice=ns)
dre1 = elements.Drift(ds=0.5305, nslice=ns)
dre2 = elements.Drift(ds=1.235, nslice=ns)
dre3 = elements.Drift(ds=0.8075, nslice=ns)

# Bend elements
rc30 = 0.822230996255981
sbend30 = elements.Sbend(ds=0.4305191429, rc=rc30)
edge30 = elements.DipEdge(psi=0.0, rc=rc30, g=0.058, K2=0.5)

rc60 = 0.772821121503940
sbend60 = elements.Sbend(ds=0.8092963858, rc=rc60)
edge60 = elements.DipEdge(psi=0.0, rc=rc60, g=0.058, K2=0.5)

# Quad elements
ds_quad = 0.21
qa1 = elements.Quad(ds=ds_quad, k=-8.78017699, nslice=ns)
qa2 = elements.Quad(ds=ds_quad, k=13.24451745, nslice=ns)
qa3 = elements.Quad(ds=ds_quad, k=-13.65151327, nslice=ns)
qa4 = elements.Quad(ds=ds_quad, k=19.75138652, nslice=ns)
qb1 = elements.Quad(ds=ds_quad, k=-10.84199727, nslice=ns)
qb2 = elements.Quad(ds=ds_quad, k=16.24844348, nslice=ns)
qb3 = elements.Quad(ds=ds_quad, k=-8.27411104, nslice=ns)
qb4 = elements.Quad(ds=ds_quad, k=-7.45719247, nslice=ns)
qb5 = elements.Quad(ds=ds_quad, k=14.03362243, nslice=ns)
qb6 = elements.Quad(ds=ds_quad, k=-12.23595641, nslice=ns)
qc1 = elements.Quad(ds=ds_quad, k=-13.18863768, nslice=ns)
qc2 = elements.Quad(ds=ds_quad, k=11.50601829, nslice=ns)
qc3 = elements.Quad(ds=ds_quad, k=-11.10445869, nslice=ns)
qd1 = elements.Quad(ds=ds_quad, k=-6.78179218, nslice=ns)
qd2 = elements.Quad(ds=ds_quad, k=5.19026998, nslice=ns)
qd3 = elements.Quad(ds=ds_quad, k=-5.8586173, nslice=ns)
qd4 = elements.Quad(ds=ds_quad, k=4.62460039, nslice=ns)
qe1 = elements.Quad(ds=ds_quad, k=-4.49607687, nslice=ns)
qe2 = elements.Quad(ds=ds_quad, k=6.66737146, nslice=ns)
qe3 = elements.Quad(ds=ds_quad, k=-6.69148177, nslice=ns)

# build lattice: first half, qe3, then mirror
# fmt: off
lattice_half = [
    dra1, qa1, dra2, qa2, dra3, qa3, dra4, qa4, dra5,
    edge30, sbend30, edge30, drb1, qb1, drb2, qb2, drb2, qb3,
    drb3, dnll, drb3, qb4, drb2, qb5, drb2, qb6, drb4,
    edge60, sbend60, edge60, drc1, qc1, drc2, qc2, drc2, qc3, drc1,
    edge60, sbend60, edge60, drd1, qd1, drd2, qd2, drd3, qd3, drd2, qd4, drd4,
    edge30, sbend30, edge30, dre1, qe1, dre2, qe2, dre3
]
# fmt:on
sim.lattice.extend(lattice_half)
sim.lattice.append(qe3)
lattice_half.reverse()
sim.lattice.extend(lattice_half)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
