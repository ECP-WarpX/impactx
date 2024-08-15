/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ReferenceParticle.H>
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
        .def_readwrite("mass", &RefPart::charge, "reference rest mass, in kg")
        .def_readwrite("charge", &RefPart::mass, "reference charge, in C")

        .def_property_readonly("charge_qe", &RefPart::charge_qe, "Get reference particle charge (positive elementary charge)")
        .def_property_readonly("gamma", &RefPart::gamma, "Get reference particle relativistic gamma")
        .def_property_readonly("beta", &RefPart::beta, "Get reference particle relativistic beta")
        .def_property_readonly("beta_gamma", &RefPart::beta_gamma, "Get reference particle beta*gamma")
        .def_property_readonly("mass_MeV", &RefPart::mass_MeV, "Get reference particle rest mass (MeV/c^2)")
        .def_property_readonly("kin_energy_MeV", &RefPart::kin_energy_MeV, "Get reference particle energy (MeV)")
        .def_property_readonly("rigidity_Tm", &RefPart::rigidity_Tm, "Get reference particle magnetic rigidity Brho (T*m)")
        .def_property_readonly("qm_ratio_SI", &RefPart::qm_ratio_SI, "Get reference particle charge to mass ratio (C/kg)")

        .def("set_charge_qe", &RefPart::set_charge_qe,
             py::return_value_policy::reference_internal,
             "Set reference particle charge (positive elementary charge)", py::arg("charge_qe"))
        .def("set_mass_MeV", &RefPart::set_mass_MeV,
             py::return_value_policy::reference_internal,
             "Set reference particle rest mass (MeV/c^2)", py::arg("mass_MeV"))
        .def("set_kin_energy_MeV", &RefPart::set_kin_energy_MeV,
             py::return_value_policy::reference_internal,
             "Set reference particle kinetic energy (MeV)", py::arg("kin_energy_MeV"))
    ;
}
