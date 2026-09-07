#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


from distribution_gsi import get_distribution, get_reference_params

from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False  # no space charge
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# set reference particle
ref = sim.beam.ref
kin_energy_MeV, bunch_charge_C, charge_qe = get_reference_params()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

# initialize the beam envelope
sim.init_envelope(ref, get_distribution(), bunch_charge_C)

# initialize accelerator lattice
ns = 10  # number of slices per ds in the element

# Drift elements
dr1 = elements.Drift(name="dr1", ds=0.6450000, nslice=ns)
dr2 = elements.Drift(name="dr2", ds=0.9700000, nslice=ns)
dr3 = elements.Drift(name="dr3", ds=6.8390117, nslice=ns)
dr4 = elements.Drift(name="dr4", ds=0.6000000, nslice=ns)
dr5 = elements.Drift(name="dr5", ds=0.7098000, nslice=ns)
dr6 = elements.Drift(name="dr6", ds=0.4998000, nslice=ns)

# Bend elements
rc = 10.0
sbend1 = elements.Sbend(name="sbend1", ds=2.6179938779914944, rc=rc)
e1 = elements.DipEdge(name="e1", psi=0.127409035395586, rc=rc, g=0.07, K2=0.5)
e2 = elements.DipEdge(name="e2", psi=0.127409035395586, rc=rc, g=0.07, K2=0.5)

# Quad elements
qs1f = elements.Quad(name="qs1f", ds=1.0400000, k=0.311872401, nslice=ns)
qs2d = elements.Quad(name="qs2d", ds=1.0400000, k=-0.496504354, nslice=ns)
qs3t = elements.Quad(name="qs3t", ds=0.4804000, k=0.62221964, nslice=ns)

# Sextupole elements
sextupole = elements.Multipole(name="sextupole", multipole=3, K_normal=0.2, K_skew=0.0)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# Lines of interest
cell = [dr1, e1, sbend1, e2, dr2, e1, sbend1, e2, dr3, qs1f, dr4, qs2d, dr5, qs3t, dr6]
chain = 11 * cell

# Construct lattice
sim.lattice.append(monitor)
sim.lattice.extend(cell)
sim.lattice.append(sextupole)  # does not affect envelope tracking
sim.lattice.extend(chain)
sim.lattice.append(monitor)

# number of turns in the ring
sim.periods = 1

# run simulation
sim.track_envelope()

# clean shutdown
sim.finalize()
