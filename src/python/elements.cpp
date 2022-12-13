/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/elements/All.H>
#include <AMReX.H>

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

    py::class_<ConstF>(me, "ConstF")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("kx"), py::arg("ky"), py::arg("kt"), py::arg("nslice") = 1,
             "A linear Constant Focusing element."
        )
        .def_property_readonly("nslice", &ConstF::nslice)
        .def_property_readonly("ds", &ConstF::ds)
    ;

    py::class_<DipEdge>(me, "DipEdge")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("psi"), py::arg("rc"), py::arg("g"), py::arg("K2"),
             "Edge focusing associated with bend entry or exit."
        )
        .def_property_readonly("nslice", &DipEdge::nslice)
        .def_property_readonly("ds", &DipEdge::ds)
    ;

    py::class_<Drift>(me, "Drift")
        .def(py::init<
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("nslice") = 1,
             "A drift."
        )
        .def_property_readonly("nslice", &Drift::nslice)
        .def_property_readonly("ds", &Drift::ds)
    ;

    py::class_<Multipole>(me, "Multipole")
        .def(py::init<
                int const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("multiple"), py::arg("K_normal"), py::arg("K_skew"),
             "A general thin multipole element."
        )
        .def_property_readonly("nslice", &Multipole::nslice)
        .def_property_readonly("ds", &Multipole::ds)
    ;

    py::class_<None>(me, "None")
        .def(py::init<>(),
             "This element does nothing."
        )
        .def_property_readonly("nslice", &None::nslice)
        .def_property_readonly("ds", &None::ds)
    ;

    py::class_<NonlinearLens>(me, "NonlinearLens")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("knll"), py::arg("cnll"),
             "Single short segment of the nonlinear magnetic insert element."
        )
        .def_property_readonly("nslice", &NonlinearLens::nslice)
        .def_property_readonly("ds", &NonlinearLens::ds)
    ;

    py::class_<Quad>(me, "Quad")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("k"), py::arg("nslice") = 1,
             "A Quadrupole magnet."
        )
        .def_property_readonly("nslice", &Quad::nslice)
        .def_property_readonly("ds", &Quad::ds)
    ;

    py::class_<RFCavity>(me, "RFCavity")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const,
                int const
             >(),
             py::arg("ds"), py::arg("escale"), py::arg("freq"),
             py::arg("phase"),
             py::arg("mapsteps"), py::arg("nslice") = 1,
             "An RF cavity (with solenoid field)."
        )
        .def_property_readonly("nslice", &RFCavity::nslice)
        .def_property_readonly("ds", &RFCavity::ds)
    ;

    py::class_<Sbend>(me, "Sbend")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("rc"), py::arg("nslice") = 1,
             "An ideal sector bend."
        )
        .def_property_readonly("nslice", &Sbend::nslice)
        .def_property_readonly("ds", &Sbend::ds)
    ;

    py::class_<ShortRF>(me, "ShortRF")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("V"), py::arg("k"),
             "A short RF cavity element at zero crossing for bunching."
        )
        .def_property_readonly("nslice", &ShortRF::nslice)
        .def_property_readonly("ds", &ShortRF::ds)
    ;
}
