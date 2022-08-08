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
        .def(py::init<>(),
             "This struct stores the reference particle attributes\n"
             "stored in ImpactXParticleContainer."
        )
        .def_readwrite("s", &RefPart::s, "integrated orbit path length, in meters")
        .def_readwrite("x", &RefPart::x, "horizontal position x, in meters")
        .def_readwrite("y", &RefPart::y, "vertical position y, in meters")
        .def_readwrite("z", &RefPart::z, "longitudinal position y, in meters")
        .def_readwrite("t", &RefPart::t, "clock time * c in meters")
        .def_readwrite("px", &RefPart::px, "momentum in x, normalized to proper velocity")
        .def_readwrite("py", &RefPart::py, "momentum in y, normalized to proper velocity")
        .def_readwrite("pz", &RefPart::pz, "momentum in z, normalized to proper velocity")
        .def_readwrite("pt", &RefPart::pt, "energy deviation, normalized by rest energy")

        .def_property_readonly("gamma", &RefPart::gamma, "Get reference particle relativistic gamma")
        .def_property_readonly("beta", &RefPart::beta, "Get reference particle relativistic beta")
        .def_property_readonly("beta_gamma", &RefPart::beta_gamma, "Get reference particle beta*gamma")
        .def("energy_MeV", &RefPart::energy_MeV, "Get reference particle energy",
             py::arg("massE_MeV"))

        .def("set_energy_MeV", &RefPart::set_energy_MeV, "Set reference particle energy",
             py::arg("energy_MeV"), py::arg("massE_MeV"))
    ;
}
