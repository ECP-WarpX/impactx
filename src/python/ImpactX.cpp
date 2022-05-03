#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <ImpactX.H>
#include <AMReX.H>

namespace py = pybind11;
using namespace impactx;

namespace {
    struct Config {};
}

void init_ImpactX(py::module& m)
{
    py::class_<ImpactX, amrex::AmrCore>(m, "ImpactX")
        .def(py::init<>())

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
