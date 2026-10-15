#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Weiqun Zhang
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
Space-charge potential of a long, uniform round beam inside a grounded pipe.

The pipe is an embedded boundary of the 3D multigrid solver. At the bunch
center, phi matches the 2D line-charge solution:

  r <= a: lambda/(4 pi eps0) (1 - r^2/a^2) + lambda/(2 pi eps0) ln(R/a)
  r >= a: lambda/(2 pi eps0) ln(R/r)
"""

import numpy as np
import pytest

from impactx import ImpactX, amr, elements

EPS0 = 8.8541878128e-12
R = 1.0e-2  # pipe radius (m)
A = 3.0e-3  # beam radius (m)
CT = 0.05  # bunch length in c*t (m)
CHARGE_C = 1.0e-9
NPART = 40000

pytestmark = pytest.mark.skipif(not amr.Config.have_eb, reason="requires ImpactX_EB=ON")


def _phi_mid_plane(**eb):
    """Track one slice; return phi in the z mid-plane, its nodes and the line charge."""
    sim = ImpactX()

    sim.n_cell = [48, 48, 32]
    sim.max_grid_size = [48]
    sim.particle_shape = 2
    sim.space_charge = "3D"
    sim.poisson_solver = "multigrid"
    sim.mlmg_relative_tolerance = 1.0e-10
    sim.prob_relative = [3.0]
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.tiny_profiler = False
    for key, value in eb.items():
        setattr(sim, "eb_" + key, value)

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(2.0e3)
    qm_eev = 1.0 / 938.27208816e6

    if amr.ParallelDescriptor.IOProcessor():
        rng = np.random.default_rng(42)
        r = A * np.sqrt(rng.random(NPART))
        theta = 2.0 * np.pi * rng.random(NPART)
        t = CT * (rng.random(NPART) - 0.5)
        zeros = np.zeros(NPART)
        sim.beam.add_n_particles(
            r * np.cos(theta),
            r * np.sin(theta),
            t,
            zeros,
            zeros,
            zeros,
            qm_eev,
            bunch_charge=CHARGE_C,
        )

    result = {}

    def save_phi(sim):
        phi = sim.phi(lev=0)
        geom = sim.Geom(lev=0)
        plane = phi[:, :, sim.n_cell[2] // 2, 0]  # MPI-collective gather
        if amr.ParallelDescriptor.IOProcessor():
            nx, ny = plane.shape
            dx = (geom.ProbHi(0) - geom.ProbLo(0)) / sim.n_cell[0]
            dy = (geom.ProbHi(1) - geom.ProbLo(1)) / sim.n_cell[1]
            result["phi"] = np.array(plane)
            result["x"] = geom.ProbLo(0) + dx * np.arange(nx)
            result["y"] = geom.ProbLo(1) + dy * np.arange(ny)
            result["dx"] = dx

    sim.hook["after_element"] = save_phi
    sim.lattice.extend([elements.Drift(name="d1", ds=1.0e-3, nslice=1)])
    sim.track_particles()

    beta = sim.beam.ref.beta
    sim.finalize()

    # line charge density in the lab frame
    result["lambda"] = CHARGE_C / (beta * CT)
    return result


def _analytic(r, lam):
    inside = lam / (4.0 * np.pi * EPS0) * (1.0 - (r / A) ** 2) + lam / (
        2.0 * np.pi * EPS0
    ) * np.log(R / A)
    outside = lam / (2.0 * np.pi * EPS0) * np.log(R / np.maximum(r, 1e-30))
    return np.where(r <= A, inside, outside)


def test_eb_circular_pipe():
    res = _phi_mid_plane(
        shape="parser",
        implicit_function=f"x^2+y^2-{R}^2",
        bounding_box_lo=[-R, -R],
        bounding_box_hi=[R, R],
    )
    if not amr.ParallelDescriptor.IOProcessor():
        return

    phi = res["phi"]
    X, Y = np.meshgrid(res["x"], res["y"], indexing="ij")
    rr = np.hypot(X, Y)
    exact = _analytic(rr, res["lambda"])
    phi0 = exact.max()

    # zero in the wall
    covered = rr > R + res["dx"]
    assert np.all(np.abs(phi[covered]) < 1.0e-8 * phi0)

    # line-charge solution between the beam edge and the wall
    gap = (rr > A + 2.0 * res["dx"]) & (rr < R - res["dx"])
    err_gap = np.max(np.abs(phi[gap] - exact[gap])) / phi0
    # and inside the beam
    core = rr < A - 2.0 * res["dx"]
    err_core = np.max(np.abs(phi[core] - exact[core])) / phi0
    print(f"max relative error: gap {err_gap:.3e}, core {err_core:.3e}")
    assert err_gap < 0.03
    assert err_core < 0.03


def test_eb_shapes_agree():
    """Named shapes match the equivalent parser or named shape."""
    ellipse = _phi_mid_plane(shape="elliptical", aperture_x=R, aperture_y=R)
    circle = _phi_mid_plane(
        shape="parser",
        implicit_function=f"x^2+y^2-{R}^2",
        bounding_box_lo=[-R, -R],
        bounding_box_hi=[R, R],
    )
    rectangle = _phi_mid_plane(shape="rectangular", aperture_x=R, aperture_y=R)
    square = _phi_mid_plane(
        shape="polygon", vertices_x=[R, -R, -R, R], vertices_y=[R, R, -R, -R]
    )
    if not amr.ParallelDescriptor.IOProcessor():
        return

    def rel_diff(a, b):
        return np.max(np.abs(a["phi"] - b["phi"])) / np.max(np.abs(b["phi"]))

    print(f"ellipse vs parser: {rel_diff(ellipse, circle):.3e}")
    print(f"rectangle vs polygon: {rel_diff(rectangle, square):.3e}")
    assert rel_diff(ellipse, circle) < 1.0e-6
    assert rel_diff(rectangle, square) < 1.0e-6
    # the square pipe has a larger aperture, so a higher potential
    assert rectangle["phi"].max() > circle["phi"].max()


if __name__ == "__main__":
    pytest.main([__file__, "-s"])
