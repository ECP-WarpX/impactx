/* Copyright 2021-2023 The ImpactX Community
 *
 * This header is used to centrally define classes that shall not violate the
 * C++ one-definition-rule (ODR) for various Python translation units.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PYIMPACTX_H
#define IMPACTX_PYIMPACTX_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include <particles/elements/All.H>

#include <list>

namespace py = pybind11;
using namespace impactx;

PYBIND11_MAKE_OPAQUE(std::list<KnownElements>)

//PYBIND11_MAKE_OPAQUE(std::vector<double>)

#endif // IMPACTX_PYIMPACTX_H
