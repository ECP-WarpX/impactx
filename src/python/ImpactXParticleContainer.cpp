/* Copyright 2022-2026 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ImpactXParticleContainer.H>
#include <diagnostics/ReducedBeamCharacteristics.H>

#include <ablastr/warn_manager/WarnManager.H>

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
        amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>
    > py_pariter_soa(m, "ImpactXParIter");
    py_pariter_soa
        .def(py::init<ParIterSoA::ContainerType&, int>(),
             py::arg("particle_container"), py::arg("level"))
        .def(py::init<ParIterSoA::ContainerType&, int, amrex::MFItInfo&>(),
             py::arg("particle_container"), py::arg("level"), py::arg("info"))
    ;

    py::class_<
        ParConstIterSoA,
        amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>
    > py_parconstiter_soa(m, "ImpactXParConstIter");
    py_parconstiter_soa
        .def(py::init<ParConstIterSoA::ContainerType&, int>(),
             py::arg("particle_container"), py::arg("level"))
        .def(py::init<ParConstIterSoA::ContainerType&, int, amrex::MFItInfo&>(),
             py::arg("particle_container"), py::arg("level"), py::arg("info"))
    ;

    py::class_<
        ImpactXParticleContainer,
        amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>
    >(m, "ImpactXParticleContainer")
        //.def(py::init<>())

        .def_property_readonly("coord_system",
            &ImpactXParticleContainer::GetCoordSystem,
            "Get the current coordinate system of particles in this container"
        )

        // simpler particle iterator loops: return types of this particle box
        // note: overwritten to return ImpactX instead of (py)AMReX iterators
        .def_property_readonly_static(
            "Iterator",
            [](py::object /* pc */){ return py::type::of<impactx::ParIterSoA>(); },
            "ImpactX iterator for particle boxes"
        )
        .def_property_readonly_static(
            "ConstIterator",
            [](py::object /* pc */){ return py::type::of<impactx::ParConstIterSoA>(); },
            "ImpactX constant iterator for particle boxes (read-only)"
        )

        .def("clear", &ImpactXParticleContainer::clear,
             py::arg("keep_mass")=false, py::arg("keep_charge")=false,
             "Empty the container and reset the reference particle"
        )
        // internal: there is a high-level Python wrapper that accepts more types
        .def("_add_n_particles",
             &ImpactXParticleContainer::AddNParticles,
             py::arg("x"), py::arg("y"), py::arg("t"),
             py::arg("px"), py::arg("py"), py::arg("pt"),
             py::arg("qm"), py::arg("bunch_charge")=py::none(), py::arg("w")=py::none(),
             py::arg("sx")=py::none(), py::arg("sy")=py::none(), py::arg("sz")=py::none(),
             "Add new particles to the container for fixed s.\n\n"
             "Either the total charge (bunch_charge) or the weight of each\n"
             "particle (w) must be provided.\n"
             "With bunch_charge, each particle represents the same number of real\n"
             "particles with the charge of the reference particle, which must be\n"
             "set up before.\n\n"
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
             ":param bunch_charge: magnitude of the total charge of the particles added in this "
             "call in C\n"
             ":param w: weight of each particle: how many real particles to represent\n"
             ":param sx: spin component in x\n"
             ":param sy: spin component in y\n"
             ":param sz: spin component in z\n"
        )
        // Getter-only property is intentional: it returns the live mutable
        // RefPart by reference.
        // A writable property (with assignment) would be ambiguous:
        // alias (Pythonic) or copy-in (safe) for the simulation-owned RefPart.
        .def_property_readonly("ref",
            [](ImpactXParticleContainer & pc) -> RefPart & {
                return pc.GetRefParticle();
            },
            "Access the reference particle."
        )
        .def("ref_particle",
            [](ImpactXParticleContainer & pc) -> RefPart & {
                py::warnings::warn(
                    "ref_particle() is deprecated. Use beam.ref instead.",
                    PyExc_DeprecationWarning,
                    2
                );
                return pc.GetRefParticle();
            },
            py::return_value_policy::reference_internal,
            "Access the reference particle.\n\n"
            "Deprecated: use ``beam.ref``."
        )
        .def("set_ref_particle",
             &ImpactXParticleContainer::SetRefParticle,
             py::arg("refpart"),
             "Set reference particle attributes."
        )
        .def("set_bucket_length",
             &ImpactXParticleContainer::SetBucketLength,
             py::arg("bucket_length"),
             "Set bucket length for particle boundary condition."
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
                 py::warnings::warn(
                    "WARNING: reduced_beam_characteristics() is deprecated. "
                    "Use beam_moments() instead.",
                    PyExc_DeprecationWarning,
                    2
                 );
                 return diagnostics::reduced_beam_characteristics(pc);
             },
             "Compute reduced beam characteristics like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters."
        )
        .def("beam_moments",
             [](ImpactXParticleContainer & pc) {
                 return pc.beam_moments();
             },
             "Calculate beam moments at current ``s`` like the position and momentum moments of the particle "
             "distribution, as well as emittance and Twiss parameters."
        )
        .def("beam_moments_history_list",
             [](ImpactXParticleContainer & pc) {
                 return pc.beam_moments_history();
             },
             "Return the history of the beam moments on every step."
        )
        .def("record_beam_moments",
             [](ImpactXParticleContainer & pc) {
                 return pc.record_beam_moments();
             },
             "Calculate & record the beam moments at current s"
        )
        .def("reset_beam_moments_history",
             [](ImpactXParticleContainer & pc) {
                 return pc.reset_beam_moments_history();
             },
             "Reset the history of the beam moments."
        )
        .def_property("store_beam_moments",
            [](ImpactXParticleContainer & pc){ return pc.store_beam_moments; },
            [](ImpactXParticleContainer & pc, bool record){ pc.store_beam_moments = record; },
            "In situ calculate and store the beam moments for every simulation step."
        )

        // TODO: cleverly pass the list of rho multifabs as references
        /*
        .def("deposit_charge",
             &ImpactXParticleContainer::DepositCharge,
             py::arg("rho"), py::arg("ref_ratio"),
             "Charge deposition"
        )
        */
    ;

    py_pariter_soa.def("pc", &ParIterSoA::pc);
    py_parconstiter_soa.def("pc", &ParConstIterSoA::pc);
}
