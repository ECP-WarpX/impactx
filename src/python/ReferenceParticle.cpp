/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ImpactXParticleContainer.H>
#include <AMReX.H>

namespace py = pybind11;
using namespace impactx;


void init_refparticle(py::module& m)
{
    py::class_<RefPart>(m, "RefPart")
        .def(py::init<>())
        .def_readwrite("x", &RefPart::x, "horizontal position x, in meters")
        .def_readwrite("y", &RefPart::y, "vertical position y, in meters")
        .def_readwrite("z", &RefPart::z, "longitudinal position y, in meters")
        .def_readwrite("t", &RefPart::t, "clock time * c in meters")
        .def_readwrite("px", &RefPart::px, "momentum in x, normalized to proper velocity")
        .def_readwrite("py", &RefPart::py, "momentum in y, normalized to proper velocity")
        .def_readwrite("pz", &RefPart::pz, "momentum in z, normalized to proper velocity")
        .def_readwrite("pt", &RefPart::pt, "energy deviation, normalized by rest energy")
    ;
}
