/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <ImpactX.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>

#if defined(AMREX_DEBUG) || defined(DEBUG)
#   include <cstdio>
#endif


namespace py = pybind11;
using namespace impactx;

namespace {
    struct Config {};
}

void init_ImpactX(py::module& m)
{
    /*
    py::bind_map<
        std::unordered_map<int, amrex::MultiFab>
    >(m, "MultiFabPerLevel");
    */

    py::class_<ImpactX>(m, "ImpactX")
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
            })
        .def("set_diags_slice_step_diagnostics",
             [](ImpactX & /* ix */, bool const slice_step_diagnostics) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("slice_step_diagnostics", slice_step_diagnostics);
             })
        .def("set_diags_file_min_digits",
             [](ImpactX & /* ix */, int const file_min_digits) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("file_min_digits", file_min_digits);
             })

        .def("init_grids", &ImpactX::initGrids)
        .def("init_beam_distribution_from_inputs", &ImpactX::initBeamDistributionFromInputs)
        .def("init_lattice_elements_from_inputs", &ImpactX::initLatticeElementsFromInputs)
        .def("evolve", &ImpactX::evolve)

        //.def_property("particle_container", &ImpactX::m_particle_container)
        //.def_readwrite("rho", &ImpactX::m_rho)
        .def_readwrite("lattice", &ImpactX::m_lattice)
        //.def_readwrite("lattice", &ImpactX::m_lattice_test)
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
