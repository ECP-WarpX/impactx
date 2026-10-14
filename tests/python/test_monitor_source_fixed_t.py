#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Test that the BeamMonitor writes and the Source reads particles at fixed t.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from impactx import (
    Config,
    CoordSystem,
    ImpactX,
    ImpactXParIter,
    RefPart,
    coordinate_transformation,
    distribution,
    elements,
    push,
)

io = pytest.importorskip("openpmd_api")

if not Config.have_openpmd:
    pytest.skip("ImpactX was compiled without openPMD support", allow_module_level=True)

# a file format that both ImpactX and the Python openPMD-api read and write
BACKEND = next(
    (
        ending
        for ending, variant in [("h5", "hdf5"), ("bp4", "adios2")]
        if Config.openpmd_backends.get(variant, False)
        and io.variants.get(variant, False)
    ),
    None,
)
if BACKEND is None:
    pytest.skip(
        "ImpactX and the Python openPMD-api share no file format",
        allow_module_level=True,
    )

# speed of light in m/s
C = 299792458.0

# tolerance in units of the rms size of the beam in each phase space coordinate
#   in single precision, the deviation from the reference momentum cancels
#   about three of the seven digits
TOL = 1.0e-10 if Config.precision == "DOUBLE" else 1.0e-3
# ... and for values that are only copied or scaled
RTOL = 1.0e-12 if Config.precision == "DOUBLE" else 1.0e-5

# the phase space coordinates at fixed s
COORDS_S = [
    "position_x",
    "position_y",
    "position_t",
    "momentum_x",
    "momentum_y",
    "momentum_t",
]

NPART = 512


def drift_to_fixed_t(x, y, t, px, py, pt, beta_gamma):
    """Particles at fixed s, at the time the reference particle is at s.

    At fixed s, t is c times the arrival time after the reference particle, px
    and py are in units of the reference momentum and pt is the energy
    deviation (gamma_ref - gamma) / beta_gamma.

    Returns the positions x, y, z relative to the reference particle and the
    momenta px, py, pz in units of the reference momentum.
    """
    gamma_ref = np.sqrt(1.0 + beta_gamma**2)

    ux = px * beta_gamma
    uy = py * beta_gamma
    gamma = gamma_ref - pt * beta_gamma
    uz = np.sqrt(gamma**2 - 1.0 - ux**2 - uy**2)

    # a particle that arrives late (t > 0) is still on its way to s
    return (
        x - ux / gamma * t,
        y - uy / gamma * t,
        -uz / gamma * t,
        px,
        py,
        uz / beta_gamma,
    )


def drift_to_fixed_s(x, y, z, ux, uy, uz, beta_gamma):
    """Particles at fixed t, as they pass the position s of the reference particle.

    At fixed t, the positions x, y, z are relative to the reference particle and
    the momenta ux, uy, uz are in units of m*c.

    Returns x, y, t, px, py, pt at fixed s.
    """
    gamma_ref = np.sqrt(1.0 + beta_gamma**2)
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)

    # c times the time until the particle reaches s: positive if it is behind
    ct = -z * gamma / uz

    return (
        x + ux / gamma * ct,
        y + uy / gamma * ct,
        ct,
        ux / beta_gamma,
        uy / beta_gamma,
        (gamma_ref - gamma) / beta_gamma,
    )


def assert_close(actual, expected, name, tol=TOL):
    """Compare a particle attribute within a fraction of its rms size."""
    scale = np.std(expected)
    if scale == 0.0:
        scale = np.max(np.abs(expected))
    np.testing.assert_allclose(
        actual, expected, rtol=0.0, atol=tol * scale, err_msg=name
    )


def new_sim(diagnostics=True):
    """A simulation without collective effects, ready for particles.

    Without diagnostics, init_grids does not move a diags directory that an
    earlier simulation wrote out of the way.
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.diagnostics = diagnostics
    sim.init_grids()

    return sim


def series_path(name, prefix="diags"):
    """The path of the openPMD series that a beam monitor wrote."""
    files = sorted(Path(prefix, "openPMD").glob(f"{name}.*"))
    assert files, f"no openPMD series found for monitor '{name}'"

    return str(files[0])


def write_monitors():
    """Track a beam to monitors at fixed s and at fixed t at the same position.

    The series at fixed s holds a second step, written after the one at fixed t.

    Returns the paths of the series at fixed s and at fixed t.
    """
    sim = new_sim()

    # a proton beam far from the speed of light, wide in angle, time and energy:
    # there, the particles at fixed s and at fixed t differ by a lot
    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(250.0)
    distr = distribution.Waterbag(
        lambdaX=1.0e-3,
        lambdaY=2.0e-3,
        lambdaT=2.0e-2,
        lambdaPx=5.0e-3,
        lambdaPy=3.0e-3,
        lambdaPt=5.0e-3,
        muxpx=-0.4,
        mutpt=0.3,
    )
    sim.add_particles(1.0e-9, distr, NPART, distribution.SpinvMF(0.4, 0.9, 0.1))

    mon_s = elements.BeamMonitor("mon_s", backend=BACKEND)
    mon_t = elements.BeamMonitor("mon_t", backend=BACKEND, coord_system="t")
    sim.lattice.extend([elements.Drift(name="d1", ds=2.0), mon_s, mon_t, mon_s])

    try:
        sim.track_particles()
    finally:
        # this closes the openPMD series, so that we can read them back below
        sim.finalize()

    return series_path("mon_s"), series_path("mon_t")


def read_step(path, step_index=-1):
    """A step of the beam in an openPMD series.

    Returns the particle attributes by name, e.g., "position_x", sorted by
    particle id, with the unitSI of each, the unitDimension of each record, the
    attributes of the beam and the time of the step.
    """
    series = io.Series(path, io.Access.read_only)
    iteration = series.iterations[sorted(series.iterations)[step_index]]
    beam = iteration.particles["beam"]

    data = {}
    unit_SI = {}
    unit_dimension = {}
    for record_name, record in beam.items():
        unit_dimension[record_name] = list(record.unit_dimension)
        for component_name, component in record.items():
            name = record_name
            if component_name != io.Record_Component.SCALAR:
                name += f"_{component_name}"
            data[name] = component.load_chunk()
            unit_SI[name] = component.unit_SI
    series.flush()

    order = np.argsort(data["id"])
    step = SimpleNamespace(
        data={name: np.array(values[order]) for name, values in data.items()},
        unit_SI=unit_SI,
        unit_dimension=unit_dimension,
        attrs={name: beam.get_attribute(name) for name in beam.attributes},
        time=iteration.time,
    )
    series.close()

    return step


def beam_particles(beam):
    """The particles in a beam: their ParticleReal attributes by name."""
    tiles = []
    for lvl in range(beam.finest_level + 1):
        for pti in ImpactXParIter(beam, level=lvl):
            tiles.append(pti.soa().to_numpy(copy=True).real)

    return {name: np.concatenate([tile[name] for tile in tiles]) for name in tiles[0]}


def assert_same_beam(actual, expected, names, tol=TOL):
    """Compare two sets of particles, which may be ordered differently."""
    assert len(actual["position_x"]) == len(expected["position_x"])

    order_actual = np.argsort(actual["position_x"])
    order_expected = np.argsort(expected["position_x"])
    for name in names:
        assert_close(
            actual[name][order_actual], expected[name][order_expected], name, tol
        )


@pytest.mark.manages_amrex
def test_monitor_fixed_t():
    """
    This tests that a beam monitor writes the particles at fixed t: drifted to
    the time of the reference particle and in the units of the openPMD standard.
    """
    path_s, path_t = write_monitors()

    at_s = read_step(path_s, 0)
    at_t = read_step(path_t)
    at_s_after = read_step(path_s, 1)

    # the output at fixed t leaves the beam that is tracked as it is
    assert at_s.data.keys() == at_s_after.data.keys()
    for name, values in at_s.data.items():
        assert np.array_equal(values, at_s_after.data[name]), name

    # the longitudinal components are named after the coordinate system
    for record_name in ["position", "momentum", "positionOffset"]:
        assert f"{record_name}_z" in at_t.data
        assert f"{record_name}_t" not in at_t.data
        assert f"{record_name}_t" in at_s.data
        assert f"{record_name}_z" not in at_s.data

    # the same particles, the same reference particle and the same beam moments
    assert len(at_t.data["id"]) == NPART
    for name in ["id", "weighting", "qm", "spin_x", "spin_y", "spin_z"]:
        assert np.array_equal(at_t.data[name], at_s.data[name]), name
    assert at_t.attrs.keys() == at_s.attrs.keys()
    for name, value in at_s.attrs.items():
        assert at_t.attrs[name] == pytest.approx(value, rel=RTOL, nan_ok=True), name
    assert np.std(at_s.data["spin_x"]) > 0.0

    # positions in meters, relative to the reference particle
    ref = at_t.attrs
    length = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    assert at_t.unit_dimension["position"] == length
    assert at_t.unit_dimension["positionOffset"] == length
    for axis in ["x", "y", "z"]:
        assert at_t.unit_SI[f"position_{axis}"] == 1.0
        assert at_t.unit_SI[f"positionOffset_{axis}"] == 1.0
        assert np.all(at_t.data[f"positionOffset_{axis}"] == ref[f"{axis}_ref"])
    assert ref["z_ref"] == pytest.approx(2.0, rel=RTOL)

    # momenta in units of the reference momentum, which unitSI converts to kg*m/s
    momentum = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    assert at_t.unit_dimension["momentum"] == momentum
    p_ref_SI = ref["beta_gamma_ref"] * ref["mass_ref"] * C
    for axis in ["x", "y", "z"]:
        assert at_t.unit_SI[f"momentum_{axis}"] == pytest.approx(p_ref_SI, rel=1.0e-12)

    # a snapshot in time, in seconds
    assert at_t.time == pytest.approx(ref["t_ref"] / C, rel=RTOL)

    # each particle drifted to the time of the reference particle
    x, y, z, px, py, pz = drift_to_fixed_t(
        *(at_s.data[name].astype(np.float64) for name in COORDS_S),
        beta_gamma=ref["beta_gamma_ref"],
    )
    assert_close(at_t.data["position_x"], x, "position_x")
    assert_close(at_t.data["position_y"], y, "position_y")
    assert_close(at_t.data["position_z"], z, "position_z")
    assert_close(at_t.data["momentum_x"], px, "momentum_x")
    assert_close(at_t.data["momentum_y"], py, "momentum_y")
    assert_close(at_t.data["momentum_z"], pz, "momentum_z")

    # ... which is not a small change for this beam
    shift_x = at_t.data["position_x"] - at_s.data["position_x"]
    assert np.std(shift_x) > 1.0e-3 * np.std(x)
    assert np.std(z) < 0.7 * np.std(at_s.data["position_t"])
    assert np.mean(at_t.data["momentum_z"]) == pytest.approx(1.0, abs=1.0e-2)


@pytest.mark.manages_amrex
def test_source_fixed_t():
    """
    This tests that a source reads the particles that a beam monitor wrote at
    fixed t and restores the beam at fixed s.
    """
    path_s, path_t = write_monitors()

    at_s = read_step(path_s)

    sim = new_sim(diagnostics=False)
    beam = sim.beam

    try:
        push(beam, elements.Source("openPMD", path_t, coord_system="t"))

        # the reference particle is restored exactly
        ref = beam.ref
        for name in ["s", "x", "y", "z", "t", "px", "py", "pz", "pt"]:
            assert getattr(ref, name) == at_s.attrs[f"{name}_ref"], name
        assert ref.kin_energy_MeV == pytest.approx(250.0, rel=RTOL)

        # ... and so are the particles
        assert_same_beam(
            beam_particles(beam),
            at_s.data,
            COORDS_S + ["weighting", "spin_x", "spin_y", "spin_z"],
        )
        assert beam.beam_moments()["charge_C"] == pytest.approx(1.0e-9, rel=RTOL)

        # a file in the other coordinate system names the option to set
        with pytest.raises(RuntimeError, match="fixed t: set coord_system to t"):
            push(beam, elements.Source("openPMD", path_t))
        with pytest.raises(RuntimeError, match="fixed s: set coord_system to s"):
            push(beam, elements.Source("openPMD", path_s, coord_system="t"))

        # the failed reads added no particles
        assert beam.total_number_of_particles() == NPART

        # both elements transform from and to a beam at fixed s
        coordinate_transformation(beam, CoordSystem.t)
        with pytest.raises(RuntimeError, match="at fixed s"):
            push(beam, elements.Source("openPMD", path_t, coord_system="t"))
        with pytest.raises(RuntimeError, match="at fixed s"):
            push(beam, elements.BeamMonitor("mon", backend=BACKEND, coord_system="t"))
        coordinate_transformation(beam, CoordSystem.s)
    finally:
        sim.finalize()


def time_based_beam(ref, z_ref):
    """A beam at fixed t as a code with time as the independent variable holds it.

    Returns positions in meters and momenta in kg*m/s in the frame of that code,
    in which the reference particle is at ref.x, ref.y, z_ref, and the number
    of real particles per macroparticle.
    """
    rng = np.random.default_rng(seed=20261013)
    npart = 400

    # in units of the reference momentum: wide in angle and momentum
    p_ref_SI = ref.beta_gamma * ref.mass * C
    return {
        "x": ref.x + 3.0e-4 + rng.normal(0.0, 1.0e-3, npart),
        "y": ref.y - 2.0e-4 + rng.normal(0.0, 2.0e-3, npart),
        "z": z_ref + rng.normal(0.0, 2.0e-2, npart),
        "px": p_ref_SI * rng.normal(0.0, 5.0e-3, npart),
        "py": p_ref_SI * rng.normal(0.0, 3.0e-3, npart),
        "pz": p_ref_SI * (1.0 + rng.normal(0.0, 5.0e-3, npart)),
        "w": rng.uniform(1.0e5, 2.0e5, npart),
    }


def write_time_based_series(path, lab, ref, layout="SI", momentum_unit=True):
    """Write a beam at fixed t as a code with time as the independent variable does.

    The layouts hold the same particles:

    - "SI": positions in meters and momenta in kg*m/s of a real particle, next to
      the charge and mass of the species, as WarpX writes them
    - "units": positions in millimeters next to a positionOffset and momenta of
      a macroparticle in units of m*c
    - "single": as "SI", in single precision
    """
    dtype = np.float32 if layout == "single" else np.float64
    npart = len(lab["x"])

    position_unit = 1.0
    momentum_unit_SI = 1.0
    offset = {"x": 0.0, "y": 0.0, "z": 0.0}
    macro_weighted = False
    if layout == "units":
        position_unit = 1.0e-3
        momentum_unit_SI = ref.mass * C
        offset = {"x": 1.0e-3, "y": -2.0e-3, "z": 0.25}
        macro_weighted = True

    series = io.Series(path, io.Access.create)
    iteration = series.iterations[100]
    iteration.time = 1.0e-9
    beam = iteration.particles["beam"]

    dataset = io.Dataset(np.dtype(dtype), [npart])
    scalar = io.Record_Component.SCALAR
    chunks = []  # openPMD reads the arrays when the series is flushed

    def store(component, values):
        chunks.append(np.ascontiguousarray(values, dtype=dtype))
        component.reset_dataset(dataset)
        component.store_chunk(chunks[-1])

    for axis in ["x", "y", "z"]:
        store(beam["position"][axis], (lab[axis] - offset[axis]) / position_unit)
        beam["position"][axis].unit_SI = position_unit

        coarse = beam["positionOffset"][axis]
        coarse.reset_dataset(dataset)
        coarse.make_constant(offset[axis] / position_unit)
        coarse.unit_SI = position_unit

        values = lab[f"p{axis}"] / momentum_unit_SI
        if macro_weighted:
            values = values * lab["w"]
        store(beam["momentum"][axis], values)
        beam["momentum"][axis].unit_SI = momentum_unit_SI

    beam["position"].unit_dimension = {io.Unit_Dimension.L: 1}
    beam["positionOffset"].unit_dimension = {io.Unit_Dimension.L: 1}
    if momentum_unit:
        beam["momentum"].unit_dimension = {
            io.Unit_Dimension.L: 1,
            io.Unit_Dimension.M: 1,
            io.Unit_Dimension.T: -1,
        }
    beam["momentum"].set_attribute("macroWeighted", np.uint32(macro_weighted))
    beam["momentum"].set_attribute("weightingPower", 1.0)

    store(beam["weighting"][scalar], lab["w"])

    for record_name, value in [("charge", ref.charge), ("mass", ref.mass)]:
        beam[record_name][scalar].reset_dataset(dataset)
        beam[record_name][scalar].make_constant(value)

    ids = np.arange(1, npart + 1, dtype=np.uint64)
    beam["id"][scalar].reset_dataset(io.Dataset(ids.dtype, [npart]))
    beam["id"][scalar].store_chunk(ids)

    series.flush()
    series.close()


@pytest.mark.manages_amrex
@pytest.mark.parametrize("layout", ["SI", "units", "single"])
def test_source_fixed_t_other_code(layout):
    """
    This tests that a source reads particles at fixed t from a code that uses
    time as the independent variable: lab frame positions and momenta in the
    units of the openPMD standard, without a reference particle in the file.
    """
    z_ref = 0.75

    sim = new_sim(diagnostics=False)
    beam = sim.beam

    try:
        # the reference particle: its species, energy and position in the file
        ref = beam.ref
        ref.set_species("proton").set_kin_energy_MeV(250.0)
        ref.x = 2.0e-4
        ref.y = -1.0e-4
        ref.z = z_ref

        lab = time_based_beam(ref, z_ref)
        path = f"time_based_{layout}.{BACKEND}"
        write_time_based_series(path, lab, ref, layout)

        push(
            beam,
            elements.Source("openPMD", path, load_ref_particle=False, coord_system="t"),
        )

        # the reference particle stays as it was configured
        assert ref.z == z_ref
        assert ref.kin_energy_MeV == pytest.approx(250.0, rel=RTOL)

        # each particle drifted to the position s of the reference particle
        x, y, t, px, py, pt = drift_to_fixed_s(
            lab["x"] - ref.x,
            lab["y"] - ref.y,
            lab["z"] - ref.z,
            lab["px"] / (ref.mass * C),
            lab["py"] / (ref.mass * C),
            lab["pz"] / (ref.mass * C),
            beta_gamma=ref.beta_gamma,
        )
        expected = dict(zip(COORDS_S, [x, y, t, px, py, pt]))
        expected["weighting"] = lab["w"]
        for axis in ["x", "y", "z"]:
            expected[f"spin_{axis}"] = np.zeros_like(x)

        # single precision in the file limits the positions and momenta
        tol = 1.0e-3 if layout == "single" else TOL
        assert_same_beam(beam_particles(beam), expected, expected.keys(), tol)

        charge_C = np.sum(lab["w"]) * ref.charge
        rtol = 1.0e-6 if layout == "single" else RTOL
        assert beam.beam_moments()["charge_C"] == pytest.approx(charge_C, rel=rtol)

        # ... which is not a small change for this beam
        assert np.std(x - (lab["x"] - ref.x)) > 1.0e-2 * np.std(x)
        assert np.std(t) > 1.4 * np.std(lab["z"])
    finally:
        sim.finalize()


@pytest.mark.manages_amrex
def test_source_fixed_t_other_code_errors():
    """
    This tests that a source names what keeps it from reading particles at
    fixed t from a code that uses time as the independent variable.
    """
    z_ref = 0.75

    sim = new_sim(diagnostics=False)
    beam = sim.beam

    try:
        ref_file = RefPart()
        ref_file.set_species("proton").set_kin_energy_MeV(250.0)
        lab = time_based_beam(ref_file, z_ref)
        path = f"time_based.{BACKEND}"
        write_time_based_series(path, lab, ref_file)

        def source(openpmd_path, **kwargs):
            return elements.Source("openPMD", openpmd_path, coord_system="t", **kwargs)

        # particles at fixed t need the momentum of the reference particle
        with pytest.raises(RuntimeError, match="no momentum"):
            push(beam, source(path, load_ref_particle=False))

        ref = beam.ref
        ref.set_species("proton").set_kin_energy_MeV(250.0)
        ref.z = z_ref

        # a file of another code holds no reference particle to load
        with pytest.raises(RuntimeError, match="load_ref_particle"):
            push(beam, source(path))

        # momenta need their unit to be told from normalized ones
        path = f"no_unit.{BACKEND}"
        write_time_based_series(path, lab, ref, momentum_unit=False)
        with pytest.raises(RuntimeError, match="momentum .* is not a momentum"):
            push(beam, source(path, load_ref_particle=False))

        # a particle that moves backward never reaches s
        lab["pz"][3] *= -1.0
        lab["pz"][7] = 0.0
        path = f"backward.{BACKEND}"
        write_time_based_series(path, lab, ref)
        with pytest.raises(RuntimeError, match="2 particles .* do not move forward"):
            push(beam, source(path, load_ref_particle=False))

        # no failed read added particles
        assert beam.total_number_of_particles() == 0
    finally:
        sim.finalize()


def run_inputs(name, inputs, diag_file_prefix):
    """Run a simulation with the beam and the lattice of an inputs file.

    Each run writes to its own directory, which keeps the output of earlier runs.
    """
    Path(name).write_text(inputs)

    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.diag_file_prefix = diag_file_prefix
    sim.init_grids()

    try:
        sim.load_inputs_file(name)
        sim.init_beam_distribution_from_inputs()
        sim.init_lattice_elements_from_inputs()
        sim.track_particles()
    finally:
        sim.finalize()


@pytest.mark.manages_amrex
def test_fixed_t_inputs_file():
    """
    This tests the coord_system option of beam_monitor and source elements in
    inputs files: one run writes the beam at fixed t and a second one reads it.
    """
    run_inputs(
        "write.in",
        f"""
beam.npart = {NPART}
beam.units = static
beam.kin_energy = 250.0
beam.charge = 1.0e-9
beam.particle = proton
beam.distribution = waterbag
beam.lambdaX = 1.0e-3
beam.lambdaY = 2.0e-3
beam.lambdaT = 2.0e-2
beam.lambdaPx = 5.0e-3
beam.lambdaPy = 3.0e-3
beam.lambdaPt = 5.0e-3

lattice.elements = drift1 mon_s mon_t

drift1.type = drift
drift1.ds = 2.0

mon_s.type = beam_monitor
mon_s.backend = {BACKEND}

mon_t.type = beam_monitor
mon_t.backend = {BACKEND}
mon_t.coord_system = t
""",
        diag_file_prefix="written",
    )
    path_t = series_path("mon_t", prefix="written")
    written = read_step(series_path("mon_s", prefix="written"))
    assert "position_z" in read_step(path_t).data

    run_inputs(
        "read.in",
        f"""
lattice.elements = source mon_s

source.type = source
source.distribution = openPMD
source.openpmd_path = "{path_t}"
source.coord_system = t

mon_s.type = beam_monitor
mon_s.backend = {BACKEND}
""",
        diag_file_prefix="read",
    )
    read = read_step(series_path("mon_s", prefix="read"))

    # the second run continues with the beam of the first one
    assert_same_beam(read.data, written.data, COORDS_S + ["weighting"])
    for name in ["s_ref", "z_ref", "t_ref", "pt_ref", "mass_ref", "charge_ref"]:
        assert read.attrs[name] == written.attrs[name], name


def test_coord_system_option():
    """
    This tests that coord_system is "s" or "t", fixed s by default, and part of
    the repr and to_dict() of both elements.
    """
    monitor = elements.BeamMonitor("monitor")
    assert monitor.coord_system == "s"
    assert monitor.to_dict()["coord_system"] == "s"

    monitor = elements.BeamMonitor("monitor", coord_system="t")
    assert monitor.coord_system == "t"

    # the dict round-trips through the constructor
    d = monitor.to_dict()
    assert d["coord_system"] == "t"
    kwargs = {k: v for k, v in d.items() if k not in ("type", "ds")}
    assert elements.BeamMonitor(**kwargs).coord_system == "t"

    with pytest.raises(RuntimeError, match="CoordSystem"):
        elements.BeamMonitor("monitor", coord_system="z")

    source = elements.Source("openPMD", "beam.h5")
    assert source.coord_system == "s"
    assert source.to_dict()["coord_system"] == "s"
    assert "coord_system=s" in repr(source)

    source = elements.Source("openPMD", "beam.h5", coord_system="t", name="source")
    assert source.coord_system == "t"
    assert source.name == "source"
    assert "coord_system=t" in repr(source)

    d = source.to_dict()
    assert d["coord_system"] == "t"
    kwargs = {k: v for k, v in d.items() if k not in ("type", "ds")}
    assert elements.Source(**kwargs).coord_system == "t"

    # the option can also be set after construction
    source.coord_system = "s"
    assert source.coord_system == "s"

    with pytest.raises(RuntimeError, match="CoordSystem"):
        elements.Source("openPMD", "beam.h5", coord_system="z")
    with pytest.raises(RuntimeError, match="CoordSystem"):
        source.coord_system = "z"
    assert source.coord_system == "s"
