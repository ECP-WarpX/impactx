/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ImpactXParticleContainer.H>
#include <AMReX.H>
#include <AMReX_ParticleContainer.H>

namespace py = pybind11;
using namespace impactx;


void init_impactxparticlecontainer(py::module& m)
{
    py::class_<
        ImpactXParticleContainer,
        amrex::ParticleContainer<0, 0, RealSoA::nattribs, IntSoA::nattribs>
    >(m, "ImpactXParticleContainer")
        //.def(py::init<>())
        .def("add_n_particles", &ImpactXParticleContainer::AddNParticles)
        .def("ref_particle",
            py::overload_cast<>(&ImpactXParticleContainer::GetRefParticle),
            py::return_value_policy::reference_internal
        )
        .def("set_ref_particle", &ImpactXParticleContainer::SetRefParticle)
        .def("set_particle_shape",
            py::overload_cast<int const>(&ImpactXParticleContainer::SetParticleShape)
        )
    ;
}
