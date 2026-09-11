/* Copyright 2022-2026 The ImpactX Community
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
#include <pybind11/functional.h>
// include <pybind11/numpy.h>  // not yet used
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/warnings.h>

#include <ImpactX.H>
#include <elements/All.H>

#include <functional>
#include <list>
#include <string>
#include <unordered_map>

namespace py = pybind11;
using namespace impactx;

PYBIND11_MAKE_OPAQUE(std::unordered_map<std::string, std::function<void(ImpactX *)> >)

#endif // IMPACTX_PYIMPACTX_H
