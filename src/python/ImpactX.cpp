/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <ImpactX.H>
#include <particles/distribution/Gaussian.H>
#include <particles/distribution/Kurth4D.H>
#include <particles/distribution/Kurth6D.H>
#include <particles/distribution/KVdist.H>
#include <particles/distribution/Semigaussian.H>
#include <particles/distribution/Waterbag.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#if defined(AMREX_DEBUG) || defined(DEBUG)
#   include <cstdio>
#endif


namespace py = pybind11;
using namespace impactx;

namespace impactx {
    struct Config {};
}

void init_ImpactX(py::module& m)
{
    py::class_<ImpactX> impactx(m, "ImpactX");
    impactx
        .def(py::init<>())

        .def("load_inputs_file",
            [](ImpactX const & /* ix */, std::string const filename) {
#if defined(AMREX_DEBUG) || defined(DEBUG)
                // note: only in debug, since this is costly for the file
                // system for highly parallel simulations with MPI
                // possible improvement:
                // - rank 0 tests file & broadcasts existence/failure
                bool inputs_file_exists = false;
                if (FILE *fp = fopen(filename.c_str(), "r")) {
                    fclose(fp);
                    inputs_file_exists = true;
                }
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(inputs_file_exists,
                    "load_inputs_file: invalid filename");
#endif

                amrex::ParmParse::addfile(filename);
            })

        .def("set_particle_shape",
            [](ImpactX & ix, int const order) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ix.m_particle_container,
                    "particle container not initialized");
                // todo: why does that not work?
                //ix.m_particle_container->SetParticleShape(order);

                amrex::ParmParse pp_ago("algo");
                pp_ago.add("particle_shape", order);
            },
            "Whether to calculate space charge effects."
        )
        .def("set_space_charge",
             [](ImpactX & /* ix */, bool const enable) {
                 amrex::ParmParse pp_algo("algo");
                 pp_algo.add("space_charge", enable);
             },
             py::arg("enable"),
             "Enable or disable space charge calculations (default: enabled)."
        )
        .def("set_diagnostics",
             [](ImpactX & /* ix */, bool const enable) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("enable", enable);
             },
             py::arg("enable"),
             "Enable or disable diagnostics generally (default: enabled).\n"
             "Disabling this is mostly used for benchmarking."
         )
        .def("set_slice_step_diagnostics",
             [](ImpactX & /* ix */, bool const enable) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("slice_step_diagnostics", enable);
             },
             py::arg("enable"),
             "Enable or disable diagnostics every slice step in elements (default: disabled).\n\n"
             "By default, diagnostics is performed at the beginning and end of the simulation.\n"
             "Enabling this flag will write diagnostics every step and slice step."
         )
        .def("set_diag_file_min_digits",
             [](ImpactX & /* ix */, int const file_min_digits) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("file_min_digits", file_min_digits);
             },
             py::arg("file_min_digits"),
             "The minimum number of digits (default: 6) used for the step\n"
             "number appended to the diagnostic file names."
         )

        .def("init_grids", &ImpactX::initGrids,
             "Initialize AMReX blocks/grids for domain decomposition & space charge mesh.\n\n"
             "This must come first, before particle beams and lattice elements are initialized."
        )
        .def("init_beam_distribution_from_inputs", &ImpactX::initBeamDistributionFromInputs)
        .def("init_lattice_elements_from_inputs", &ImpactX::initLatticeElementsFromInputs)
        .def("add_particles", &ImpactX::add_particles,
             py::arg("qm"), py::arg("bunch_charge"),
             py::arg("distr"), py::arg("npart"),
             "Generate and add n particles to the particle container.\n\n"
             "Will also resize the geometry based on the updated particle\n"
             "distribution's extent and then redistribute particles in according\n"
             "AMReX grid boxes."
        )
        .def("evolve", &ImpactX::evolve,
             "Run the main simulation loop for a number of steps."
        )
        .def("particle_container",
             [](ImpactX & ix) -> ImpactXParticleContainer & {
                return *ix.m_particle_container;
             },
             py::return_value_policy::reference_internal,
             "Access the beam particle container."
        )
        //.def_readwrite("rho", &ImpactX::m_rho)
        .def_readwrite("lattice",
            &ImpactX::m_lattice,
            "Access the accelerator element lattice."
        )
    ;

    py::class_<Config>(m, "Config")
//        .def_property_readonly_static(
//            "impactx_version",
//            [](py::object) { return Version(); },
//            "ImpactX version")
        .def_property_readonly_static(
            "have_mpi",
            [](py::object){
#ifdef AMREX_USE_MPI
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "have_gpu",
            [](py::object){
#ifdef AMREX_USE_GPU
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "have_omp",
            [](py::object){
#ifdef AMREX_USE_OMP
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "gpu_backend",
            [](py::object){
#ifdef AMREX_USE_CUDA
                return "CUDA";
#elif defined(AMREX_USE_HIP)
                return "HIP";
#elif defined(AMREX_USE_DPCPP)
                return "SYCL";
#else
                return py::none();
#endif
            })
        ;
}
