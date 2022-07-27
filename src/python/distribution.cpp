/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/distribution/All.H>

#include <AMReX.H>
#include <AMReX_REAL.H>

#include <variant>

namespace py = pybind11;
using namespace impactx;


void init_distribution(py::module& m)
{
    py::module_ md = m.def_submodule(
        "distribution",
        "Particle beam distributions in ImpactX"
    );

    py::class_<distribution::Gaussian>(md, "Gaussian")
        .def(py::init<
                amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0
        );

    py::class_<distribution::Kurth4D>(md, "Kurth4D")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
         >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0
        );

    py::class_<distribution::Kurth6D>(md, "Kurth6D")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0
        );

    py::class_<distribution::KVdist>(md, "KVdist")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0
        );

    py::class_<distribution::Semigaussian>(md, "Semigaussian")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0
        );

    py::class_<distribution::Waterbag>(md, "Waterbag")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0
        );
}
