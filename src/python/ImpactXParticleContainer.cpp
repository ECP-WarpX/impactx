/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ImpactXParticleContainer.H>
#include <particles/diagnostics/ReducedBeamCharacteristics.H>

#include <AMReX.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParticleContainer.H>

#include <algorithm>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace impactx;


void init_impactxparticlecontainer(py::module& m)
{
    py::enum_<CoordSystem>(m, "CoordSystem")
        .value("s", CoordSystem::s)
        .value("t", CoordSystem::t)
        .export_values();

    py::class_<
        ParIterSoA,
        amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs>
    >(m, "ImpactXParIter")
        .def(py::init<ParIterSoA::ContainerType&, int>(),
             py::arg("particle_container"), py::arg("level"))
        .def(py::init<ParIterSoA::ContainerType&, int, amrex::MFItInfo&>(),
             py::arg("particle_container"), py::arg("level"), py::arg("info"))
    ;

    py::class_<
        ParConstIterSoA,
        amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs>
    >(m, "ImpactXParConstIter")
        .def(py::init<ParConstIterSoA::ContainerType&, int>(),
             py::arg("particle_container"), py::arg("level"))
        .def(py::init<ParConstIterSoA::ContainerType&, int, amrex::MFItInfo&>(),
             py::arg("particle_container"), py::arg("level"), py::arg("info"))
    ;

    py::class_<
        ImpactXParticleContainer,
        amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs>
    >(m, "ImpactXParticleContainer")
        //.def(py::init<>())

        .def_property_readonly_static("RealSoA",
            [](py::object /* pc */){ return py::type::of<RealSoA>(); },
            "RealSoA attribute name labels"
        )

        .def_property_readonly("coord_system",
            &ImpactXParticleContainer::GetCoordSystem,
            "Get the current coordinate system of particles in this container"
        )

        // simpler particle iterator loops: return types of this particle box
        // note: overwritten to return ImpactX instead of (py)AMReX iterators
        .def_property_readonly_static(
            "iterator",
            [](py::object /* pc */){ return py::type::of<impactx::ParIterSoA>(); },
            "ImpactX iterator for particle boxes"
        )
        .def_property_readonly_static(
            "const_iterator",
            [](py::object /* pc */){ return py::type::of<impactx::ParConstIterSoA>(); },
            "ImpactX constant iterator for particle boxes (read-only)"
        )

        .def("add_n_particles",
             &ImpactXParticleContainer::AddNParticles,
             py::arg("x"), py::arg("y"), py::arg("t"),
             py::arg("px"), py::arg("py"), py::arg("pt"),
             py::arg("qm"), py::arg("bchchg"),
             "Add new particles to the container for fixed s.\n\n"
             "Note: This can only be used *after* the initialization (grids) have\n"
             "      been created, meaning after the call to ImpactX.init_grids\n"
             "      has been made in the ImpactX class.\n\n"
             ":param x: positions in x\n"
             ":param y: positions in y\n"
             ":param t: positions as time-of-flight in c*t\n"
             ":param px: momentum in x\n"
             ":param py: momentum in y\n"
             ":param pt: momentum in t\n"
             ":param qm: charge over mass in 1/eV\n"
             ":param bchchg: total charge within a bunch in C"
        )
        .def("ref_particle",
            py::overload_cast<>(&ImpactXParticleContainer::GetRefParticle),
            py::return_value_policy::reference_internal,
            "Access the reference particle."
        )
        .def("set_ref_particle",
             &ImpactXParticleContainer::SetRefParticle,
             py::arg("refpart"),
             "Set reference particle attributes."
        )
        .def("min_and_max_positions",
             &ImpactXParticleContainer::MinAndMaxPositions,
             "Compute the min and max of the particle position in each dimension.\n\n"
             ":return: x_min, y_min, z_min, x_max, y_max, z_max"
        )
        .def("mean_and_std_positions",
             &ImpactXParticleContainer::MeanAndStdPositions,
             "Compute the mean and std of the particle position in each dimension.\n\n"
             ":return: x_mean, x_std, y_mean, y_std, z_mean, z_std"
        )
        .def("reduced_beam_characteristics",
             [](ImpactXParticleContainer & pc) {
                 return diagnostics::reduced_beam_characteristics(pc);
             },
             "Compute reduced beam characteristics like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters."
        )

        .def("redistribute",
             &ImpactXParticleContainer::Redistribute,
             "Redistribute particles in the current mesh in x, y, z"
        )
        // TODO: cleverly pass the list of rho multifabs as references
        /*
        .def("deposit_charge",
             &ImpactXParticleContainer::DepositCharge,
             py::arg("rho"), py::arg("ref_ratio"),
             "Charge deposition"
        )
        */

        .def_property_readonly("RealSoA_names", &ImpactXParticleContainer::RealSoA_names,
              "Get the name of each ParticleReal SoA component")
        .def_property_readonly("intSoA_names", &ImpactXParticleContainer::intSoA_names,
               "Get the name of each int SoA component")
    ;

    m.def("get_RealSoA_names", &get_RealSoA_names,
          py::arg("num_real_comps"),
          "Get the name of each ParticleReal SoA component\n\nnum_real_comps: pass number of compile-time + runtime arrays");
    m.def("get_intSoA_names", &get_intSoA_names,
          py::arg("num_int_comps"),
          "Get the name of each int SoA component\n\nnum_int_comps: pass number of compile-time + runtime arrays");
}
