/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/elements/All.H>
#include <AMReX.H>

#include <vector>

namespace py = pybind11;
using namespace impactx;


void init_elements(py::module& m)
{
    py::module_ me = m.def_submodule(
        "elements",
        "Accelerator lattice elements in ImpactX"
    );

    using KnownElementsList = std::list<KnownElements>;
    py::class_<KnownElementsList> kel(me, "KnownElementsList");
    kel
        .def(py::init<>())
        .def(py::init<KnownElements>())
        .def(py::init([](py::list l){
            auto v = new KnownElementsList;
            for (auto const & handle : l)
                v->push_back(handle.cast<KnownElements>());
            return v;
        }))

        .def("append", [](KnownElementsList &v, KnownElements el) { v.emplace_back(el); },
             "Add a single element to the list.")

        .def("extend",
             [](KnownElementsList &v, KnownElementsList l) {
                 for (auto const & el : l)
                    v.push_back(el);
                 return v;
             },
             "Add a list of elements to the list.")
        .def("extend",
             [](KnownElementsList &v, py::list l) {
                 for (auto const & handle : l)
                 {
                    auto el = handle.cast<KnownElements>();
                    v.push_back(el);
                 }
                 return v;
             },
             "Add a list of elements to the list."
         )

        .def("clear", &KnownElementsList::clear,
             "Clear the list to become empty.")
        .def("pop_back", &KnownElementsList::pop_back,
             "Return and remove the last element of the list.")
        .def("__len__", [](const KnownElementsList &v) { return v.size(); },
             "The length of the list.")
        .def("__iter__", [](KnownElementsList &v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>()) /* Keep list alive while iterator is used */
    ;

    // mixin classes

    py::class_<elements::Thick>(me, "Thick")
        .def(py::init<
                 amrex::ParticleReal const,
                 amrex::ParticleReal const
             >(),
             py::arg("ds"), py::arg("nslice") = 1,
             "Mixin class for lattice elements with finite length."
        )
        .def_property_readonly("nslice", &elements::Thick::nslice)
        .def_property_readonly("ds", &elements::Thick::ds)
    ;

    py::class_<elements::Thin>(me, "Thin")
        .def(py::init<>(),
             "Mixin class for lattice elements with zero length."
        )
        .def_property_readonly("nslice", &elements::Thin::nslice)
        .def_property_readonly("ds", &elements::Thin::ds)
    ;

    // beam optics below

    py::class_<ConstF, elements::Thick>(me, "ConstF")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("kx"), py::arg("ky"), py::arg("kt"), py::arg("nslice") = 1,
             "A linear Constant Focusing element."
        )
    ;

    py::class_<DipEdge, elements::Thin>(me, "DipEdge")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("psi"), py::arg("rc"), py::arg("g"), py::arg("K2"),
             "Edge focusing associated with bend entry or exit."
        )
    ;

    py::class_<Drift, elements::Thick>(me, "Drift")
        .def(py::init<
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("nslice") = 1,
             "A drift."
        )
    ;

    py::class_<Multipole, elements::Thin>(me, "Multipole")
        .def(py::init<
                int const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("multiple"), py::arg("K_normal"), py::arg("K_skew"),
             "A general thin multipole element."
        )
    ;

    py::class_<None, elements::Thin>(me, "None")
        .def(py::init<>(),
             "This element does nothing."
        )
    ;

    py::class_<NonlinearLens, elements::Thin>(me, "NonlinearLens")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("knll"), py::arg("cnll"),
             "Single short segment of the nonlinear magnetic insert element."
        )
    ;

    py::class_<Programmable>(me, "Programmable")
        .def(py::init<>(),
             "A programmable beam optics element."
        )
        .def_property("nslice",
            [](Programmable & p) { return p.nslice(); },
            [](Programmable & p, int nslice) { p.m_nslice = nslice; }
        )
        .def_property("ds",
              [](Programmable & p) { return p.ds(); },
              [](Programmable & p, amrex::ParticleReal ds) { p.m_ds = ds; }
        )
        .def_property("beam_particles",
              [](Programmable & p) { return p.m_beam_particles; },
              [](Programmable & p,
                 std::function<void(ImpactXParticleContainer::iterator *, RefPart &)> new_hook
              ) { p.m_beam_particles = new_hook; }
        )
        .def_property("ref_particle",
              [](Programmable & p) { return p.m_ref_particle; },
              [](Programmable & p,
                 std::function<void(RefPart &)> new_hook
              ) { p.m_ref_particle = new_hook; }
        )
    ;

    py::class_<Quad, elements::Thick>(me, "Quad")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("k"), py::arg("nslice") = 1,
             "A Quadrupole magnet."
        )
    ;

    py::class_<RFCavity, elements::Thick>(me, "RFCavity")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                std::vector<amrex::ParticleReal>,
                std::vector<amrex::ParticleReal>,
                int const,
                int const
             >(),
             py::arg("ds"), py::arg("escale"), py::arg("freq"),
             py::arg("phase"), py::arg("cos_coefficients"), py::arg("sin_coefficients"),
             py::arg("mapsteps") = 1, py::arg("nslice") = 1,
             "An RF cavity."
        )
    ;

    py::class_<Sbend, elements::Thick>(me, "Sbend")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("rc"), py::arg("nslice") = 1,
             "An ideal sector bend."
        )
    ;

    py::class_<ShortRF, elements::Thin>(me, "ShortRF")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("V"), py::arg("k"),
             "A short RF cavity element at zero crossing for bunching."
        )
    ;

    py::class_<SoftSolenoid, elements::Thick>(me, "SoftSolenoid")
        .def(py::init<
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 std::vector<amrex::ParticleReal>,
                 std::vector<amrex::ParticleReal>,
                 int,
                 int
             >(),
             py::arg("ds"), py::arg("bscale"),
             py::arg("cos_coefficients"), py::arg("sin_coefficients"),
             py::arg("mapsteps") = 1, py::arg("nslice") = 1,
             "A soft-edge solenoid."
        )
    ;

    py::class_<Sol, elements::Thick>(me, "Sol")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("ks"), py::arg("nslice") = 1,
             "An ideal hard-edge Solenoid magnet."
        )
    ;
}
