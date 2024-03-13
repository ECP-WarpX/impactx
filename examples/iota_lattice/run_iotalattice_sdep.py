#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# diagnostics: IOTA nonlinear lens invariants calculation
sim.set_diag_iota_invariants(
    alpha=1.376381920471173, beta=1.892632003628881, tn=0.4, cn=0.01
)

# init particle beam
energy_MeV = 2.5
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=1.865379469388e-003,
    lambdaY=2.0192133150418e-003,
    lambdaT=1.0e-3,
    lambdaPx=1.402566720991e-003,
    lambdaPy=9.57593913381e-004,
    lambdaPt=0.0,
    muxpx=0.482260919078473,
    muypy=0.762127656873158,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

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

# Special (elliptic) nonlinear element:

# set basic parameters of the nonlinear element
lens_length = 1.8
num_lenses = 18
tune_advance = 0.3
c_parameter = 0.01
t_strength = 0.4
ds = lens_length / num_lenses

# build up the nonlinear lens in segments
ds_half = ds / 2.0
dr = elements.Drift(ds=ds_half, nslice=1)
nll = []
for j in range(0, num_lenses):
    s = -lens_length / 2.0 + ds_half + j * ds
    beta_star = lens_length / 2.0 * 1.0 / math.tan(math.pi * tune_advance)
    beta = beta_star * (
        1.0 + (2.0 * s * math.tan(math.pi * tune_advance) / lens_length) ** 2
    )
    knll_s = t_strength * c_parameter**2 * ds / beta
    cnll_s = c_parameter * math.sqrt(beta)
    nllens = elements.NonlinearLens(knll=knll_s, cnll=cnll_s)
    segment = [dr, nllens, dr]
    nll.extend(segment)

lattice_before_nll = [
    dra1,
    qa1,
    dra2,
    qa2,
    dra3,
    qa3,
    dra4,
    qa4,
    dra5,
    edge30,
    sbend30,
    edge30,
    drb1,
    qb1,
    drb2,
    qb2,
    drb2,
    qb3,
    drb3,
]

lattice_after_nll = [
    drb3,
    qb4,
    drb2,
    qb5,
    drb2,
    qb6,
    drb4,
    edge60,
    sbend60,
    edge60,
    drc1,
    qc1,
    drc2,
    qc2,
    drc2,
    qc3,
    drc1,
    edge60,
    sbend60,
    edge60,
    drd1,
    qd1,
    drd2,
    qd2,
    drd3,
    qd3,
    drd2,
    qd4,
    drd4,
    edge30,
    sbend30,
    edge30,
    dre1,
    qe1,
    dre2,
    qe2,
    dre3,
]

# build lattice: first half, qe3, then mirror
# modified to place nll as the first element
# fmt:on
sim.lattice.append(monitor)
sim.lattice.extend(nll)
sim.lattice.extend(lattice_after_nll)
sim.lattice.append(qe3)
lattice_after_nll.reverse()
sim.lattice.extend(lattice_after_nll)
sim.lattice.append(dnll)
lattice_before_nll.reverse()
sim.lattice.extend(lattice_before_nll)
lattice_before_nll.reverse()
sim.lattice.extend(lattice_before_nll)
sim.lattice.append(monitor)

# number of turns in the ring
sim.periods = 5

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
