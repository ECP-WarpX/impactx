#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Tests for ``unroll_impactx_lattice`` of the optional ``synmadx`` MAD-X parser.
"""

import pytest

import impactx
from impactx import elements

synmadx = pytest.importorskip("impactx.synmadx")


def test_unroll_keeps_unnamed_elements_unnamed():
    """An element without a name is recreated without a name, not named "None"."""
    lattice = [
        elements.Drift(ds=1.0),
        elements.Quad(ds=0.5, k=1.0, name="q1"),
        elements.Sbend(ds=0.5, rc=10.0),
        elements.DipEdge(psi=0.1, rc=10.0, g=0.0, K2=0.0),
    ]

    source = synmadx.unroll_impactx_lattice(lattice)
    rebuilt = eval(source, {"impactx": impactx})

    assert [type(e).__name__ for e in rebuilt] == [type(e).__name__ for e in lattice]
    assert [e.name for e in rebuilt] == [None, "q1", None, None]
    assert [e.has_name for e in rebuilt] == [False, True, False, False]
