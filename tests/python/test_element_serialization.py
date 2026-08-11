#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Test that all elements can be serialized via to_dicts() and reconstructed via from_dicts().
"""

import re

import numpy as np
import pytest

import amrex.space3d as amr
from impactx import Config, elements

# amrex Real SmallMatrix precision suffix matching the build
_REAL = "float" if Config.precision == "SINGLE" else "double"


def get_constructor_params(element_class):
    """
    Extract constructor parameter names from pybind11 docstring.

    Returns a set of parameter names, or None if parsing fails.
    """
    doc = element_class.__init__.__doc__
    if not doc:
        return None

    # Find the signature line: __init__(self: ..., param1: type, ...) -> None
    match = re.search(r"__init__\(self:[^,]+,\s*(.+?)\)\s*->", doc)
    if not match:
        match = re.search(r"__init__\(self:[^,]+,\s*(.+?)\)\s*$", doc, re.MULTILINE)
    if not match:
        match = re.search(r"__init__\(self:[^)]+\)", doc)
        if match:
            return set()
        return None

    params_str = match.group(1) if match.lastindex else ""
    if not params_str.strip():
        return set()

    # Parse individual parameters, handling nested types like dict[str, int]
    params = []
    depth = 0
    current = ""
    for char in params_str:
        if char in "[(":
            depth += 1
        elif char in "])":
            depth -= 1
        elif char == "," and depth == 0:
            if current.strip():
                params.append(current.strip())
            current = ""
            continue
        current += char
    if current.strip():
        params.append(current.strip())

    param_names = set()
    for p in params:
        name = p.split(":")[0].strip()
        if name:
            param_names.add(name)

    return param_names


def dicts_equal(d1: dict, d2: dict, rtol: float = 1e-12) -> bool:
    """
    Compare two element dicts for equality.

    Handles floating point comparison with relative tolerance and
    compares AMReX SmallMatrix objects element-wise.
    """
    if set(d1.keys()) != set(d2.keys()):
        return False

    for key in d1:
        v1, v2 = d1[key], d2[key]

        if isinstance(v1, float) and isinstance(v2, float):
            if v1 == 0.0 and v2 == 0.0:
                continue
            if v1 == 0.0 or v2 == 0.0:
                return False
            if abs(v1 - v2) / max(abs(v1), abs(v2)) > rtol:
                return False
        elif isinstance(v1, list) and isinstance(v2, list):
            if len(v1) != len(v2):
                return False
            for a, b in zip(v1, v2):
                if isinstance(a, float) and isinstance(b, float):
                    if a == 0.0 and b == 0.0:
                        continue
                    if a == 0.0 or b == 0.0:
                        return False
                    if abs(a - b) / max(abs(a), abs(b)) > rtol:
                        return False
                elif a != b:
                    return False
        elif hasattr(v1, "to_numpy") and hasattr(v2, "to_numpy"):
            arr1 = v1.to_numpy()
            arr2 = v2.to_numpy()
            if not np.allclose(arr1, arr2, rtol=rtol, atol=0):
                return False
        elif v1 != v2:
            return False

    return True


# Elements that cannot be roundtripped
SKIP_ELEMENTS = {
    "Source",  # Requires a valid openPMD file path
    "Empty",  # BUG: to_dict returns type='None' instead of 'Empty'
}


@pytest.fixture
def all_elements():
    """Create a KnownElementsList with one instance of each element type."""
    lattice = elements.KnownElementsList()

    lattice.append(
        elements.Aperture(
            aperture_x=0.01,
            aperture_y=0.02,
            repeat_x=0.0,
            repeat_y=0.0,
            shape="rectangular",
            action="transmit",
            dx=0.001,
            dy=0.002,
            rotation=0.1,
            name="test_aperture",
        )
    )

    monitor = elements.BeamMonitor(
        name="test_monitor",
        backend="h5",
        encoding="g",
        period_sample_intervals=2,
    )
    lattice.append(monitor)

    monitor_nopart = elements.BeamMonitor(
        name="monitor_nopart",
        backend="h5",
        encoding="f",
        period_sample_intervals=3,
        particles=False,
        beam_moments=True,
    )
    lattice.append(monitor_nopart)

    lattice.append(
        elements.Buncher(
            V=0.5,
            k=10.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.03,
            name="test_buncher",
        )
    )

    lattice.append(
        elements.CFbend(
            ds=0.5,
            rc=2.0,
            k=0.5,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            nslice=3,
            name="test_cfbend",
        )
    )

    lattice.append(
        elements.ChrAcc(
            ds=1.0,
            ez=1e6,
            bz=0.5,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            nslice=2,
            name="test_chracc",
        )
    )

    lattice.append(
        elements.ChrDrift(
            ds=1.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            nslice=2,
            name="test_chrdrift",
        )
    )

    lattice.append(
        elements.ChrPlasmaLens(
            ds=0.2,
            k=100.0,
            unit=0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.01,
            aperture_y=0.01,
            nslice=2,
            name="test_chrplasmalens",
        )
    )

    lattice.append(
        elements.ChrQuad(
            ds=0.3,
            k=2.0,
            unit=0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            nslice=2,
            name="test_chrquad",
        )
    )

    lattice.append(
        elements.ConstF(
            ds=0.5,
            kx=1.0,
            ky=1.5,
            kt=0.5,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            nslice=2,
            name="test_constf",
        )
    )

    lattice.append(
        elements.DipEdge(
            psi=0.1,
            rc=5.0,
            g=0.01,
            K2=0.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.03,
            name="test_dipedge",
        )
    )

    lattice.append(
        elements.Drift(
            ds=1.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            nslice=2,
            name="test_drift",
        )
    )

    lattice.append(
        elements.ExactCFbend(
            ds=0.5,
            k_normal=[1.0, 0.5],
            k_skew=[0.0, 0.1],
            unit=0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            int_order=4,
            mapsteps=5,
            nslice=2,
            name="test_exactcfbend",
        )
    )

    lattice.append(
        elements.ExactDrift(
            ds=1.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            nslice=2,
            name="test_exactdrift",
        )
    )

    lattice.append(
        elements.ExactMultipole(
            ds=0.5,
            k_normal=[1.0, 0.5],
            k_skew=[0.0, 0.1],
            unit=0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            int_order=4,
            mapsteps=5,
            nslice=2,
            name="test_exactmultipole",
        )
    )

    lattice.append(
        elements.ExactQuad(
            ds=0.3,
            k=2.0,
            unit=0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            nslice=2,
            mapsteps=10,
            int_order=4,
            name="test_exactquad",
        )
    )

    lattice.append(
        elements.ExactSbend(
            ds=0.5,
            phi=15.0,
            B=0.5,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            nslice=2,
            name="test_exactsbend",
        )
    )

    lattice.append(
        elements.Kicker(
            xkick=0.001,
            ykick=0.002,
            unit="dimensionless",
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            name="test_kicker",
        )
    )

    R = getattr(amr, f"SmallMatrix_6x6_F_SI1_{_REAL}").identity()
    lattice.append(
        elements.LinearMap(
            R=R,
            ds=0.5,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.03,
            name="test_linearmap",
        )
    )

    lattice.append(elements.Marker(name="test_marker"))

    lattice.append(
        elements.Multipole(
            multipole=2,
            K_normal=0.1,
            K_skew=0.05,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.03,
            name="test_multipole",
        )
    )

    lattice.append(
        elements.NonlinearLens(
            knll=0.5,
            cnll=0.01,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            name="test_nonlinearlens",
        )
    )

    lattice.append(
        elements.PlaneXYRot(
            angle=45.0, dx=0.001, dy=0.002, rotation=0.05, name="test_planexyrot"
        )
    )

    vertices_x = [0.01, 0.01, -0.01, -0.01, 0.01]
    vertices_y = [0.01, -0.01, -0.01, 0.01, 0.01]
    lattice.append(
        elements.PolygonAperture(
            vertices_x=vertices_x,
            vertices_y=vertices_y,
            min_radius2=1e-6,
            repeat_x=0.0,
            repeat_y=0.0,
            action="transmit",
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            name="test_polygonaperture",
        )
    )

    lattice.append(elements.Programmable(ds=0.5, nslice=2, name="test_programmable"))

    lattice.append(elements.PRot(phi_in=0.1, phi_out=-0.1, name="test_prot"))

    lattice.append(
        elements.Quad(
            ds=0.3,
            k=2.5,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            nslice=2,
            name="test_quad",
        )
    )

    lattice.append(
        elements.QuadEdge(
            k=2.0,
            unit=0,
            flag="entry",
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.03,
            name="test_quadedge",
        )
    )

    lattice.append(
        elements.RFCavity(
            ds=0.5,
            escale=1e6,
            freq=1.3e9,
            phase=-90.0,
            cos_coefficients=[1.0],
            sin_coefficients=[0.0],
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            nslice=2,
            mapsteps=10,
            name="test_rfcavity",
        )
    )

    lattice.append(
        elements.Sbend(
            ds=0.5,
            rc=5.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            nslice=2,
            name="test_sbend",
        )
    )

    lattice.append(
        elements.ShortRF(
            V=0.5,
            freq=1.3e9,
            phase=-90.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            name="test_shortrf",
        )
    )

    lattice.append(
        elements.Sol(
            ds=0.5,
            ks=1.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            nslice=2,
            name="test_sol",
        )
    )

    lattice.append(
        elements.SoftQuadrupole(
            ds=0.5,
            gscale=1.0,
            cos_coefficients=[1.0, 0.5],
            sin_coefficients=[0.0, 0.1],
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            mapsteps=5,
            nslice=2,
            name="test_softquadrupole",
        )
    )

    lattice.append(
        elements.SoftSolenoid(
            ds=0.5,
            bscale=1.0,
            cos_coefficients=[1.0, 0.5],
            sin_coefficients=[0.0, 0.1],
            unit=0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.02,
            mapsteps=5,
            nslice=2,
            name="test_softsolenoid",
        )
    )

    v = getattr(amr, f"SmallMatrix_3x1_F_SI1_{_REAL}")()
    v.set_val(0.1)
    A = getattr(amr, f"SmallMatrix_3x6_F_SI1_{_REAL}")()
    A.set_val(0.01)
    lattice.append(
        elements.SpinMap(
            v=v,
            A=A,
            ds=0.5,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            name="test_spinmap",
        )
    )

    lattice.append(
        elements.TaperedPL(
            k=100.0,
            taper=10.0,
            unit=0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.02,
            aperture_y=0.03,
            name="test_taperedpl",
        )
    )

    lattice.append(
        elements.ThinDipole(
            theta=0.05,
            rc=10.0,
            dx=0.001,
            dy=0.002,
            rotation=0.05,
            aperture_x=0.03,
            aperture_y=0.04,
            name="test_thindipole",
        )
    )

    yield lattice, monitor
    monitor.finalize()


def test_to_dicts_from_dicts_roundtrip(all_elements):
    """
    Test that lattice serialization via to_dicts()/from_dicts() roundtrips correctly.
    """
    lattice, _ = all_elements

    # Serialize
    dicts = lattice.to_dicts()
    assert len(dicts) == len(lattice)

    # Verify each dict has required 'type' key
    for d in dicts:
        assert "type" in d

    # Reconstruct
    lattice2 = elements.KnownElementsList()
    lattice2.from_dicts(dicts)
    assert len(lattice2) == len(lattice)

    # Serialize again and compare
    dicts2 = lattice2.to_dicts()

    for i, (d1, d2) in enumerate(zip(dicts, dicts2)):
        assert dicts_equal(d1, d2), (
            f"Dict mismatch at index {i} ({d1.get('type')}):\n"
            f"  Original: {d1}\n"
            f"  Reconstructed: {d2}"
        )


def test_element_dict_has_constructor_params(all_elements):
    """
    Test that to_dict() returns all constructor parameters for each element.
    """
    lattice, _ = all_elements

    for element in lattice:
        element_type = type(element).__name__

        if element_type in SKIP_ELEMENTS:
            continue

        d = element.to_dict()
        constructor_params = get_constructor_params(type(element))

        if constructor_params is None:
            continue

        dict_keys = set(d.keys())
        dict_keys.discard("type")
        if d.get("ds") == 0.0 and "ds" not in constructor_params:
            dict_keys.discard("ds")

        missing = constructor_params - dict_keys
        extra = dict_keys - constructor_params

        assert not missing, f"{element_type}: to_dict() missing params: {missing}"
        assert not extra, f"{element_type}: to_dict() has extra keys: {extra}"


def test_mixin_properties_are_settable(all_elements):
    """
    Test that the read/write properties inherited from the mixin classes can be
    assigned after construction, and that to_dict() sees the new value.
    """
    lattice, _ = all_elements

    # property name -> value to assign
    settable = {
        "aperture_x": 0.05,
        "aperture_y": 0.06,
        "dx": 0.007,
        "dy": 0.008,
        "rotation": 0.09,
    }
    # apertures are stored inverted, so a set/get pair costs two divisions
    rtol = 1e-5 if Config.precision == "SINGLE" else 1e-12

    covered = set()
    for element in lattice:
        element_type = type(element).__name__

        if element_type in SKIP_ELEMENTS:
            continue

        for prop, value in settable.items():
            if not hasattr(element, prop):
                continue

            setattr(element, prop, value)
            assert getattr(element, prop) == pytest.approx(value, rel=rtol), (
                f"{element_type}.{prop} did not keep the assigned value"
            )

            d = element.to_dict()
            assert prop in d, f"{element_type}: to_dict() is missing {prop}"
            assert d[prop] == pytest.approx(value, rel=rtol), (
                f"{element_type}: to_dict()[{prop!r}] does not reflect the assignment"
            )

            covered.add(prop)

    assert covered == set(settable), (
        f"no element exercised: {sorted(set(settable) - covered)}"
    )


def test_aperture_disabled_by_nonpositive():
    """
    Test that Aperture.aperture_x/aperture_y follow the PipeAperture convention:
    a half-aperture of zero or less disables the boundary in that plane and
    reads back as zero, both from the constructor and from assignment.
    """
    aperture = elements.Aperture(aperture_x=0.01, aperture_y=0.02)

    # a valid assignment is kept and is visible in to_dict()
    aperture.aperture_x = 0.03
    aperture.aperture_y = 0.04
    assert aperture.aperture_x == pytest.approx(0.03)
    assert aperture.aperture_y == pytest.approx(0.04)
    assert aperture.to_dict()["aperture_x"] == pytest.approx(0.03)
    assert aperture.to_dict()["aperture_y"] == pytest.approx(0.04)

    # zero or less disables the boundary and reads back as an exact zero,
    # so that a to_dict() roundtrip reproduces the disabled element
    for disabled in [0.0, -0.01]:
        aperture.aperture_x = disabled
        aperture.aperture_y = disabled
        assert aperture.aperture_x == 0.0
        assert aperture.aperture_y == 0.0
        assert aperture.to_dict()["aperture_x"] == 0.0
        assert aperture.to_dict()["aperture_y"] == 0.0

        # the constructor uses the same convention
        from_ctor = elements.Aperture(aperture_x=disabled, aperture_y=disabled)
        assert from_ctor.aperture_x == 0.0
        assert from_ctor.aperture_y == 0.0


def test_individual_element_roundtrip(all_elements):
    """
    Test each element type individually for roundtrip consistency.
    """
    lattice, _ = all_elements

    for element in lattice:
        element_type = type(element).__name__

        if element_type in SKIP_ELEMENTS:
            continue

        # work-around for .to_dict() returning radians
        # where contructor expects degrees
        buggy_rad_angle_types = [
            "ExactSbend",
            "PlaneXYRot",
            "PRot",
            "ThinDipole",
        ]
        extra_kwargs = {}
        if element_type in buggy_rad_angle_types:
            extra_kwargs = {"in_degrees": True}

        # Roundtrip single element
        d1 = element.to_dict(**extra_kwargs)
        temp_lattice = elements.KnownElementsList()
        temp_lattice.from_dicts([d1])
        element2 = temp_lattice[0]
        d2 = element2.to_dict(**extra_kwargs)

        assert dicts_equal(d1, d2), (
            f"Roundtrip failed for {element_type}:\n"
            f"  Original: {d1}\n"
            f"  Reconstructed: {d2}"
        )


def test_to_py_roundtrip(all_elements):
    """
    Test that to_py() generates executable code that recreates the lattice.
    """
    lattice, _ = all_elements

    # Generate Python code
    code = lattice.to_py()

    # Verify code structure
    assert "from impactx import elements" in code
    assert "def get_lattice():" in code
    assert "return lattice" in code

    # Execute the generated code
    namespace = {}
    exec(code, namespace)

    # Get the reconstructed lattice
    get_lattice = namespace["get_lattice"]
    lattice2 = get_lattice()

    assert len(lattice2) == len(lattice)

    # Compare all elements
    dicts1 = lattice.to_dicts()
    dicts2 = lattice2.to_dicts()

    for i, (d1, d2) in enumerate(zip(dicts1, dicts2)):
        assert dicts_equal(d1, d2), (
            f"to_py() roundtrip failed at index {i} ({d1.get('type')}):\n"
            f"  Original: {d1}\n"
            f"  Reconstructed: {d2}"
        )


def test_copy_covers_all_element_types(all_elements):
    """Every element type must survive ``copy()`` with its configuration intact.

    ``replace_each()`` copies the template it is given once per selected position, so a
    type that cannot be copied would break that operation for any lattice containing one.
    """
    lattice, _ = all_elements

    failures = []
    for element in lattice:
        type_name = type(element).__name__
        if type_name in SKIP_ELEMENTS:
            continue
        try:
            duplicate = element.copy()
        except Exception as e:  # noqa: BLE001 - report every offender at once
            failures.append(f"{type_name}: {type(e).__name__}: {e}")
            continue
        assert type(duplicate) is type(element), type_name
        assert dicts_equal(duplicate.to_dict(), element.to_dict()), (
            f"copy of {type_name} does not match the original:\n"
            f"  Original: {element.to_dict()}\n"
            f"  Copy:     {duplicate.to_dict()}"
        )

    assert not failures, "elements that cannot be cloned:\n  " + "\n  ".join(failures)


def test_lattice_rebuild_covers_all_element_types(all_elements):
    """The same coverage through the public API that depends on cloning."""
    lattice, _ = all_elements

    n_before = len(lattice)
    # delete a single element: every *other* element has to be cloned
    lattice.select(kind="Marker").delete()
    assert len(lattice) < n_before
