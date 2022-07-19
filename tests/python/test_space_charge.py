# -*- coding: utf-8 -*-

from impactx import ImpactX, elements
import matplotlib.pyplot as plt


def tst_space_charge():
    """
    This tests an equivalent to main.cpp in C++
    """
    impactX = ImpactX()

    impactX.load_inputs_file("examples/kurth/input_kurth.in")

    impactX.init_grids()
    impactX.init_beam_distribution_from_inputs()
    impactX.init_lattice_elements_from_inputs()

    impactX.evolve(num_steps=1)

    # rho = impactX...
