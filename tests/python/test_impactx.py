# -*- coding: utf-8 -*-

import impactx


def test_impactx():
    """
    This tests an equivalent to main.cpp in C++
    """
    impactX = impactx.ImpactX()

    impactX.init_grids()

    # TODO: not yet working to add runtime files; work in AMReX needed
    # needs https://github.com/AMReX-Codes/amrex/pull/2842
    impactX.load_inputs_file("examples/fodo/input_fodo.in")

    impactX.init_beam_distribution_from_inputs()
    impactX.init_lattice_elements_from_inputs()

    impactX.evolve(num_steps=1)
