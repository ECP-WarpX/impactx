# -*- coding: utf-8 -*-

from impactx import ImpactX, elements


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
    impactX.set_diags_slice_step_diagnostics(True)
    impactX.init_grids()

    # init particle beam
    # TODO

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

    # TODO: enable once particle beam is loaded
    #impactX.evolve()
