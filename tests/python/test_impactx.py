# -*- coding: utf-8 -*-

from impactx import ImpactX, RefPart, distribution, elements


def test_impactx_fodo_file():
    """
    This tests an equivalent to main.cpp in C++
    """
    impactX = ImpactX()

    impactX.load_inputs_file("examples/fodo/input_fodo.in")

    impactX.init_grids()
    impactX.init_beam_distribution_from_inputs()
    impactX.init_lattice_elements_from_inputs()

    impactX.evolve()


def test_impactx_nofile():
    """
    This tests using ImpactX without an inputs file
    """
    impactX = ImpactX()

    impactX.set_particle_shape(2)
    impactX.set_slice_step_diagnostics(True)
    impactX.init_grids()

    # init particle beam
    energy_MeV = 2.0e3
    charge_C = 0.0
    mass_MeV = 0.510998950
    qm_qeeV = -1.0/0.510998950e6
    npart = 10000

    distr = distribution.Waterbag(
        sigmaX = 3.9984884770e-5,
        sigmaY = 3.9984884770e-5,
        sigmaT = 1.0e-3,
        sigmaPx = 2.6623538760e-5,
        sigmaPy = 2.6623538760e-5,
        sigmaPt = 2.0e-3,
        muxpx = -0.846574929020762,
        muypy = 0.846574929020762,
        mutpt = 0.0)
    impactX.add_particles(qm_qeeV, charge_C, distr, npart)

    # init reference particle
    refPart = RefPart()
    # make the next two lines a helper function?
    refPart.pt = -energy_MeV / mass_MeV - 1.0
    refPart.pz = (refPart.pt**2 - 1.0)**0.5
    impactX.particle_container().set_ref_particle(refPart)

    assert(impactX.particle_container().TotalNumberOfParticles() == npart)

    # init accelerator lattice
    fodo = [
        elements.Drift(0.25),
        elements.Quad(1.0, 1.0),
        elements.Drift(0.5),
        elements.Quad(1.0, -1.0),
        elements.Drift(0.25)
    ]
    #  assign a fodo segment
    #impactX.lattice = fodo

    #  add 4 more FODO segments
    for i in range(4):
        impactX.lattice.extend(fodo)

    # add 2 more drifts
    for i in range(4):
        impactX.lattice.append(elements.Drift(0.25))

    print(impactX.lattice)
    print(len(impactX.lattice))
    assert(len(impactX.lattice) > 5)

    impactX.evolve()
