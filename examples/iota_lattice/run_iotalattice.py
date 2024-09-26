#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# init particle beam
kin_energy_MeV = 2.5
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=1.588960728035e-3,
    lambdaY=2.496625268437e-3,
    lambdaT=1.0e-3,
    lambdaPx=2.8320397837724e-3,
    lambdaPy=1.802433091137e-3,
    lambdaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# init accelerator lattice
ns = 10  # number of slices per ds in the element

# Drift elements
dra1 = elements.Drift(name="dra1", ds=0.9125, nslice=ns)
dra2 = elements.Drift(name="dra2", ds=0.135, nslice=ns)
dra3 = elements.Drift(name="dra3", ds=0.725, nslice=ns)
dra4 = elements.Drift(name="dra4", ds=0.145, nslice=ns)
dra5 = elements.Drift(name="dra5", ds=0.3405, nslice=ns)
drb1 = elements.Drift(name="drb1", ds=0.3205, nslice=ns)
drb2 = elements.Drift(name="drb2", ds=0.14, nslice=ns)
drb3 = elements.Drift(name="drb3", ds=0.1525, nslice=ns)
drb4 = elements.Drift(name="drb4", ds=0.31437095, nslice=ns)
drc1 = elements.Drift(name="drc1", ds=0.42437095, nslice=ns)
drc2 = elements.Drift(name="drc2", ds=0.355, nslice=ns)
dnll = elements.Drift(name="dnll", ds=1.8, nslice=ns)
drd1 = elements.Drift(name="drd1", ds=0.62437095, nslice=ns)
drd2 = elements.Drift(name="drd2", ds=0.42, nslice=ns)
drd3 = elements.Drift(name="drd3", ds=1.625, nslice=ns)
drd4 = elements.Drift(name="drd4", ds=0.6305, nslice=ns)
dre1 = elements.Drift(name="dre1", ds=0.5305, nslice=ns)
dre2 = elements.Drift(name="dre2", ds=1.235, nslice=ns)
dre3 = elements.Drift(name="dre3", ds=0.8075, nslice=ns)

# Bend elements
rc30 = 0.822230996255981
sbend30 = elements.Sbend(name="sbend30", ds=0.4305191429, rc=rc30)
edge30 = elements.DipEdge(name="edge30", psi=0.0, rc=rc30, g=0.058, K2=0.5)

rc60 = 0.772821121503940
sbend60 = elements.Sbend(name="sbend60", ds=0.8092963858, rc=rc60)
edge60 = elements.DipEdge(name="edge60", psi=0.0, rc=rc60, g=0.058, K2=0.5)

# Quad elements
ds_quad = 0.21
qa1 = elements.Quad(name="qa1", ds=ds_quad, k=-8.78017699, nslice=ns)
qa2 = elements.Quad(name="qa2", ds=ds_quad, k=13.24451745, nslice=ns)
qa3 = elements.Quad(name="qa3", ds=ds_quad, k=-13.65151327, nslice=ns)
qa4 = elements.Quad(name="qa4", ds=ds_quad, k=19.75138652, nslice=ns)
qb1 = elements.Quad(name="qb1", ds=ds_quad, k=-10.84199727, nslice=ns)
qb2 = elements.Quad(name="qb2", ds=ds_quad, k=16.24844348, nslice=ns)
qb3 = elements.Quad(name="qb3", ds=ds_quad, k=-8.27411104, nslice=ns)
qb4 = elements.Quad(name="qb4", ds=ds_quad, k=-7.45719247, nslice=ns)
qb5 = elements.Quad(name="qb5", ds=ds_quad, k=14.03362243, nslice=ns)
qb6 = elements.Quad(name="qb6", ds=ds_quad, k=-12.23595641, nslice=ns)
qc1 = elements.Quad(name="qc1", ds=ds_quad, k=-13.18863768, nslice=ns)
qc2 = elements.Quad(name="qc2", ds=ds_quad, k=11.50601829, nslice=ns)
qc3 = elements.Quad(name="qc3", ds=ds_quad, k=-11.10445869, nslice=ns)
qd1 = elements.Quad(name="qd1", ds=ds_quad, k=-6.78179218, nslice=ns)
qd2 = elements.Quad(name="qd2", ds=ds_quad, k=5.19026998, nslice=ns)
qd3 = elements.Quad(name="qd3", ds=ds_quad, k=-5.8586173, nslice=ns)
qd4 = elements.Quad(name="qd4", ds=ds_quad, k=4.62460039, nslice=ns)
qe1 = elements.Quad(name="qe1", ds=ds_quad, k=-4.49607687, nslice=ns)
qe2 = elements.Quad(name="qe2", ds=ds_quad, k=6.66737146, nslice=ns)
qe3 = elements.Quad(name="qe3", ds=ds_quad, k=-6.69148177, nslice=ns)

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
sim.lattice.append(monitor)
sim.lattice.extend(lattice_half)
sim.lattice.append(qe3)
lattice_half.reverse()
sim.lattice.extend(lattice_half)
sim.lattice.append(monitor)

# number of turns in the ring
sim.periods = 5

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
