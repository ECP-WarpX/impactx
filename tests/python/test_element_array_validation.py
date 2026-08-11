#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
"""Array-valued element parameters are validated before anything reads them."""

import pytest

from impactx import elements

# elements whose field evaluation reads the first coefficient, so they need at least one
COEFFICIENT_ELEMENTS = {
    "ExactCFbend": (
        dict(ds=0.1, k_normal=[1.0], k_skew=[0.0]),
        ("k_normal", "k_skew"),
    ),
    "RFCavity": (
        dict(
            ds=0.1,
            escale=1.0,
            freq=1.0e9,
            phase=0.0,
            cos_coefficients=[2.0],
            sin_coefficients=[0.0],
        ),
        ("cos_coefficients", "sin_coefficients"),
    ),
    "SoftQuadrupole": (
        dict(ds=0.1, gscale=1.0, cos_coefficients=[2.0], sin_coefficients=[0.0]),
        ("cos_coefficients", "sin_coefficients"),
    ),
    "SoftSolenoid": (
        dict(ds=0.1, bscale=1.0, cos_coefficients=[2.0], sin_coefficients=[0.0]),
        ("cos_coefficients", "sin_coefficients"),
    ),
}


def build(name):
    kwargs, _ = COEFFICIENT_ELEMENTS[name]
    return getattr(elements, name)(**dict(kwargs))


@pytest.mark.parametrize("name", sorted(COEFFICIENT_ELEMENTS))
def test_emptying_the_coefficients_is_rejected(name):
    """An element that reads a coefficient may not be left with none."""

    element = build(name)

    with pytest.raises(ValueError):
        element.set_coefficients([], [])


@pytest.mark.parametrize("name", sorted(COEFFICIENT_ELEMENTS))
def test_constructing_without_coefficients_is_rejected(name):
    kwargs, (first, second) = COEFFICIENT_ELEMENTS[name]
    kwargs = dict(kwargs)
    kwargs[first] = []
    kwargs[second] = []

    with pytest.raises(ValueError):
        getattr(elements, name)(**kwargs)


@pytest.mark.parametrize("name", sorted(COEFFICIENT_ELEMENTS))
def test_mismatched_lengths_are_rejected(name):
    element = build(name)

    with pytest.raises(ValueError):
        element.set_coefficients([1.0, 2.0], [1.0])


def test_multipole_transfer_map_survives_a_single_coefficient():
    """A pure dipole has no quadrupole component and supplies no second coefficient."""

    from impactx import RefPart

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    element = elements.ExactMultipole(ds=0.1, k_normal=[1.0], k_skew=[0.0])
    element.transfer_map(ref)  # must not read past the end

    element.set_coefficients([2.0], [0.0])
    element.transfer_map(ref)


class TestPolygonVertices:
    OUTLINE_X = [-1.0, 1.0, 1.0, -1.0, -1.0]
    OUTLINE_Y = [-1.0, -1.0, 1.0, 1.0, -1.0]

    def polygon(self):
        return elements.PolygonAperture(
            vertices_x=list(self.OUTLINE_X), vertices_y=list(self.OUTLINE_Y)
        )

    def test_an_empty_partner_is_rejected_not_indexed(self):
        """The lengths are checked before the outline is read positionally."""

        polygon = self.polygon()

        with pytest.raises(ValueError, match="same length"):
            polygon.set_vertices(list(self.OUTLINE_X), [])

    def test_mismatched_lengths_say_so(self):
        """The message names the real problem rather than the closed-outline rule."""

        polygon = self.polygon()

        with pytest.raises(ValueError, match="same length"):
            polygon.set_vertices([0.0, 1.0, 0.0], [0.0, 1.0])

    def test_an_open_outline_is_rejected(self):
        polygon = self.polygon()

        with pytest.raises(ValueError, match="first and last vertex"):
            polygon.set_vertices([0.0, 1.0, 2.0], [0.0, 1.0, 2.0])

    def test_a_rejected_update_leaves_the_vertices_alone(self):
        polygon = self.polygon()

        with pytest.raises(ValueError):
            polygon.set_vertices(list(self.OUTLINE_X), [])

        assert list(polygon.vertices_x) == self.OUTLINE_X
        assert list(polygon.vertices_y) == self.OUTLINE_Y
