#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
    py::class_<ImpactX, amrex::AmrCore>(m, "ImpactX")
        .def(py::init<>())

        .def("load_inputs_file",
            [](ImpactX const & /* ix */, std::string const filename) {
#if defined(AMREX_DEBUG) || defined(DEBUG)
                // note: only in debug, since this is costly for the file
                // system for highly parallel simulations with MPI
                // possible improvement:
                // - rank 0 tests file & broadcasts existance/failure
                bool inputs_file_exists = false;
                if (FILE *fp = fopen(filename.c_str(), "r")) {
                    fclose(fp);
                    inputs_file_exists = true;
                }
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(inputs_file_exists,
                    "load_inputs_file: invalid filename");
#endif

                // TODO: needs https://github.com/AMReX-Codes/amrex/pull/2842
                amrex::ParmParse::addfile(filename);
            })
        .def("init_grids", &ImpactX::initGrids)
        .def("init_beam_distribution_from_inputs", &ImpactX::initBeamDistributionFromInputs)
        .def("init_lattice_elements_from_inputs", &ImpactX::initLatticeElementsFromInputs)
        .def("evolve", &ImpactX::evolve, py::arg("num_steps"))
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
