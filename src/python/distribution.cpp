/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <initialization/InitDistribution.H>
#include <particles/distribution/All.H>

#include <AMReX.H>
#include <AMReX_REAL.H>

#include <variant>

namespace py = pybind11;
using namespace impactx;

namespace
{
    /** Register impactx::generate_add_particles for each distribution
     *
     * @tparam I counter through all known distributions in impactx::distribution::KnownDistributions
     * @param md the module to register the function in
     */
    template< std::size_t I = 0 >
    void register_generate_add_particles(py::module& md)
    {
        using V = impactx::distribution::KnownDistributions;
        using T = std::variant_alternative_t<I, V>;

        md.def("generate_add_particles", &generate_add_particles<T>);

        if constexpr (I < std::variant_size_v<V> - 1)
            register_generate_add_particles<I + 1>(md);
    }
}

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
             py::arg("muxpx"), py::arg("muypy"), py::arg("mutpt")
        );

    py::class_<distribution::Kurth4D>(md, "Kurth4D")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
         >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx"), py::arg("muypy"), py::arg("mutpt")
        );

    py::class_<distribution::Kurth6D>(md, "Kurth6D")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx"), py::arg("muypy"), py::arg("mutpt")
        );

    py::class_<distribution::KVdist>(md, "KVdist")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx"), py::arg("muypy"), py::arg("mutpt")
        );

    py::class_<distribution::Semigaussian>(md, "Semigaussian")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx"), py::arg("muypy"), py::arg("mutpt")
        );

    py::class_<distribution::Waterbag>(md, "Waterbag")
        .def(py::init<
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const,
                 amrex::ParticleReal const, amrex::ParticleReal const, amrex::ParticleReal const
             >(),
             py::arg("sigmaX"), py::arg("sigmaY"), py::arg("sigmaT"),
             py::arg("sigmaPx"), py::arg("sigmaPy"), py::arg("sigmaPt"),
             py::arg("muxpx"), py::arg("muypy"), py::arg("mutpt")
        );

    register_generate_add_particles(md);
}
