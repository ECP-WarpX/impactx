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
            for (auto &handle : l)
                v->push_back(handle.cast<KnownElements>());
            return v;
        }))

        .def("append", [](KnownElementsList &v, KnownElements el) { v.emplace_back(el); })

        .def("extend", [](KnownElementsList &v, KnownElementsList l) {
            for (auto &el : l)
                v.push_back(el);
            return v;
        })
        .def("extend", [](KnownElementsList &v, py::list l) {
            for (auto &handle : l)
            {
                auto el = handle.cast<KnownElements>();
                v.push_back(el);
            }
            return v;
        })

        .def("clear", &KnownElementsList::clear)
        .def("pop_back", &KnownElementsList::pop_back)
        .def("__len__", [](const KnownElementsList &v) { return v.size(); })
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
             py::arg("ds"), py::arg("kx"), py::arg("ky"), py::arg("kt"), py::arg("nslice") = 1
        )
        .def_property_readonly("nslice", &ConstF::nslice)
    ;

    py::class_<DipEdge>(me, "DipEdge")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("psi"), py::arg("rc"), py::arg("g"), py::arg("K2")
        )
        .def_property_readonly("nslice", &DipEdge::nslice)
    ;

    py::class_<Drift>(me, "Drift")
        .def(py::init<
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("nslice") = 1
        )
        .def_property_readonly("nslice", &Drift::nslice)
    ;

    py::class_<Multipole>(me, "Multipole")
        .def(py::init<
                int const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("multiple"), py::arg("K_normal"), py::arg("K_skew")
        )
        .def_property_readonly("nslice", &Multipole::nslice)
    ;

    py::class_<NonlinearLens>(me, "NonlinearLens")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("knll"), py::arg("cnll")
        )
        .def_property_readonly("nslice", &NonlinearLens::nslice)
    ;

    py::class_<Sbend>(me, "Sbend")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("rc"), py::arg("nslice") = 1
        )
        .def_property_readonly("nslice", &Sbend::nslice)
    ;

    py::class_<ShortRF>(me, "ShortRF")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("V"), py::arg("k")
        )
        .def_property_readonly("nslice", &ShortRF::nslice)
    ;

    py::class_<Quad>(me, "Quad")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("k"), py::arg("nslice") = 1
        )
        .def_property_readonly("nslice", &Quad::nslice)
    ;
}
