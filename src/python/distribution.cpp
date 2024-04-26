/* Copyright 2021-2023 The ImpactX Community
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
    py::module_ const md = m.def_submodule(
        "distribution",
        "Particle beam distributions in ImpactX"
    );

    py::class_<distribution::Gaussian>(md, "Gaussian")
        .def(py::init<
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambdaX"), py::arg("lambdaY"), py::arg("lambdaT"),
             py::arg("lambdaPx"), py::arg("lambdaPy"), py::arg("lambdaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0,
             "A 6D Gaussian distribution"
        );

    py::class_<distribution::Kurth4D>(md, "Kurth4D")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
         >(),
             py::arg("lambdaX"), py::arg("lambdaY"), py::arg("lambdaT"),
             py::arg("lambdaPx"), py::arg("lambdaPy"), py::arg("lambdaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0,
             "A 4D Kurth distribution transversely + a uniform distribution\n"
             "in t + a Gaussian distribution in pt"
        );

    py::class_<distribution::Kurth6D>(md, "Kurth6D")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambdaX"), py::arg("lambdaY"), py::arg("lambdaT"),
             py::arg("lambdaPx"), py::arg("lambdaPy"), py::arg("lambdaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0,
             "A 6D Kurth distribution\n\n"
             "R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)\n"
             "C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)"
        );

    py::class_<distribution::KVdist>(md, "KVdist")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambdaX"), py::arg("lambdaY"), py::arg("lambdaT"),
             py::arg("lambdaPx"), py::arg("lambdaPy"), py::arg("lambdaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0,
             "A K-V distribution transversely + a uniform distribution\n"
             "in t + a Gaussian distribution in pt"
        );

    py::class_<distribution::None>(md, "None")
        .def(py::init<>(),
             "Sets all values to zero."
        );

    py::class_<distribution::Semigaussian>(md, "Semigaussian")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambdaX"), py::arg("lambdaY"), py::arg("lambdaT"),
             py::arg("lambdaPx"), py::arg("lambdaPy"), py::arg("lambdaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0,
             "A 6D Semi-Gaussian distribution (uniform in position, Gaussian in momentum)."
        );

    py::class_<distribution::Thermal>(md, "Thermal")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("k"), py::arg("kT"), py::arg("kT_halo"),
             py::arg("normalize"), py::arg("normalize_halo"),
             py::arg("halo")=0.0,
             "A stationary thermal or bithermal distribution\n\n"
             "R. D. Ryne, J. Qiang, and A. Adelmann, in Proc. EPAC2004, pp. 1942-1944 (2004)"
        );

    py::class_<distribution::Triangle>(md, "Triangle")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambdaX"), py::arg("lambdaY"), py::arg("lambdaT"),
             py::arg("lambdaPx"), py::arg("lambdaPy"), py::arg("lambdaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0,
             "A triangle distribution for laser-plasma acceleration related applications.\n\n"
             "A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).\n"
             "The transverse distribution is a 4D waterbag."
        );

    py::class_<distribution::Waterbag>(md, "Waterbag")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambdaX"), py::arg("lambdaY"), py::arg("lambdaT"),
             py::arg("lambdaPx"), py::arg("lambdaPy"), py::arg("lambdaPt"),
             py::arg("muxpx")=0.0, py::arg("muypy")=0.0, py::arg("mutpt")=0.0,
             "A 6D Waterbag distribution"
        );
}
