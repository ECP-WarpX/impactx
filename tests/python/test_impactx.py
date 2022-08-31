# -*- coding: utf-8 -*-

from impactx import ImpactX, RefPart, distribution, elements


def test_impactx_fodo_file():
    """
    This tests an equivalent to main.cpp in C++
    """
    sim = ImpactX()

    sim.load_inputs_file("examples/fodo/input_fodo.in")

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()

    sim.evolve()


def test_impactx_nofile():
    """
    This tests using ImpactX without an inputs file
    """
    sim = ImpactX()

    sim.set_particle_shape(2)
    sim.set_slice_step_diagnostics(True)
    sim.init_grids()

    # init particle beam
    energy_MeV = 2.0e3
    bunch_charge_C = 1.0e-9
    npart = 10000

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

    assert sim.particle_container().TotalNumberOfParticles() == npart

    # init accelerator lattice
    fodo = [
        elements.Drift(0.25),
        elements.Quad(1.0, 1.0),
        elements.Drift(0.5),
        elements.Quad(1.0, -1.0),
        elements.Drift(0.25),
    ]
    #  assign a fodo segment
    # sim.lattice = fodo

    #  add 4 more FODO segments
    for i in range(4):
        sim.lattice.extend(fodo)

    # add 2 more drifts
    for i in range(4):
        sim.lattice.append(elements.Drift(0.25))

    print(sim.lattice)
    print(len(sim.lattice))
    assert len(sim.lattice) > 5

    sim.evolve()
