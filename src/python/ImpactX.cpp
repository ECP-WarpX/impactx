/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <ImpactX.H>
#include <particles/diagnostics/ReducedBeamCharacteristics.H>
#include <particles/distribution/Gaussian.H>
#include <particles/distribution/Kurth4D.H>
#include <particles/distribution/Kurth6D.H>
#include <particles/distribution/KVdist.H>
#include <particles/distribution/Semigaussian.H>
#include <particles/distribution/Waterbag.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#if defined(AMREX_DEBUG) || defined(DEBUG)
#   include <cstdio>
#endif
#include <string>


namespace py = pybind11;
using namespace impactx;

namespace impactx {
    struct Config {};
}

namespace detail
{
    /** Helper Function for Property Getters
     *
     * This queries an amrex::ParmParse entry. This throws a
     * std::runtime_error if the entry is not found.
     *
     * This handles most common throw exception logic in ImpactX instead of
     * going over library boundaries via amrex::Abort().
     *
     * @tparam T type of the amrex::ParmParse entry
     * @param prefix the prefix, e.g., "impactx" or "amr"
     * @param name the actual key of the entry, e.g., "particle_shape"
     * @return the queried value (or throws if not found)
     */
    template< typename T>
    auto get_or_throw (std::string const prefix, std::string const name)
    {
        T value;
        bool const has_name = amrex::ParmParse(prefix).query(name.c_str(), value);

        if (!has_name)
            throw std::runtime_error(prefix + "." + name + " is not set yet");
        return value;
    }
}

void init_ImpactX (py::module& m)
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
                    "load_inputs_file: file does not exist: " + filename);
#endif

                amrex::ParmParse::addfile(filename);
            })

        .def_property("n_cell",
            [](ImpactX & /* ix */) {
                std::vector<int> n_cell;
                amrex::ParmParse pp_amr("amr");
                pp_amr.getarr("n_cell", n_cell);
                return n_cell;
            },
            [](ImpactX & /* ix */, std::array<int, AMREX_SPACEDIM> /* n_cell */) {
                throw std::runtime_error("setting n_cell is not yet implemented");
                /*
                amrex::ParmParse pp_amr("amr");
                amrex::Vector<int> const n_cell_v(n_cell.begin(), n_cell.end());
                pp_amr.addarr("n_cell", n_cell_v);

                int const max_level = ix.maxLevel();
                for (int lev=0; lev<=max_level; lev++) {
                    ix.ClearLevel(lev);
                    // TODO: more...
                }
                if (amrex::ParallelDescriptor::IOProcessor())
                    ix.printGridSummary(amrex::OutStream(), 0, max_level);
                */
            },
            "The number of grid points along each direction on the coarsest level."
        )

        .def_property("domain",
            [](ImpactX & /* ix */) {
                amrex::ParmParse pp_geometry("geometry");
                amrex::Vector<amrex::Real> prob_lo;
                amrex::Vector<amrex::Real> prob_hi;
                pp_geometry.getarr("prob_lo", prob_lo);
                pp_geometry.getarr("prob_hi", prob_hi);
                amrex::RealBox rb(prob_lo.data(), prob_hi.data());
                return rb;
            },
            [](ImpactX & ix, amrex::RealBox rb) {
                amrex::ParmParse pp_geometry("geometry");
                amrex::RealVect const prob_lo_rv(rb.lo());
                amrex::Vector<amrex::Real> const prob_lo_v({rb.lo()[0], rb.lo()[1], rb.lo()[2]});
                amrex::Vector<amrex::Real> const prob_hi_v({rb.hi()[0], rb.hi()[1], rb.hi()[2]});
                pp_geometry.addarr("prob_lo", prob_lo_v);
                pp_geometry.addarr("prob_hi", prob_hi_v);

                pp_geometry.add("dynamic_size", false);

                ix.ResizeMesh();
            },
            "The physical extent of the full simulation domain, relative to the reference particle position, in meters."
        )

        .def_property("prob_relative",
              [](ImpactX & /* ix */) {
                  return detail::get_or_throw<amrex::Real>("geometry", "prob_relative");
              },
              [](ImpactX & /* ix */, amrex::Real frac) {
                  amrex::ParmParse pp_geometry("geometry");
                  pp_geometry.add("prob_relative", frac);
              },
              "The field mesh is expanded beyond the physical extent of particles by this factor."
        )

        .def_property("dynamic_size",
              [](ImpactX & /* ix */) {
                  amrex::ParmParse pp_geometry("geometry");
                  bool dynamic_size;
                  pp_geometry.get("dynamic_size", dynamic_size);
                  return dynamic_size;
              },
              [](ImpactX & /* ix */, bool dynamic_size) {
                  amrex::ParmParse pp_geometry("geometry");
                  pp_geometry.add("dynamic_size", dynamic_size);
              },
              "Use dynamic (``true``) resizing of the field mesh or static sizing (``false``)."
        )

        .def_property("particle_shape",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<int>("algo", "particle_shape");
            },
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
        .def_property("space_charge",
             [](ImpactX & /* ix */) {
                 return detail::get_or_throw<bool>("algo", "space_charge");
             },
             [](ImpactX & /* ix */, bool const enable) {
                 amrex::ParmParse pp_algo("algo");
                 pp_algo.add("space_charge", enable);
             },
             "Enable or disable space charge calculations (default: enabled)."
        )
        .def_property("diagnostics",
             [](ImpactX & /* ix */) {
                 return detail::get_or_throw<bool>("diag", "enable");
             },
             [](ImpactX & /* ix */, bool const enable) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("enable", enable);
             },
             "Enable or disable diagnostics generally (default: enabled).\n"
             "Disabling this is mostly used for benchmarking."
         )
        .def_property("slice_step_diagnostics",
             [](ImpactX & /* ix */) {
                 return detail::get_or_throw<bool>("diag", "slice_step_diagnostics");
             },
             [](ImpactX & /* ix */, bool const enable) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("slice_step_diagnostics", enable);
             },
             "Enable or disable diagnostics every slice step in elements (default: disabled).\n\n"
             "By default, diagnostics is performed at the beginning and end of the simulation.\n"
             "Enabling this flag will write diagnostics every step and slice step."
         )
        .def_property("diag_file_min_digits",
             [](ImpactX & /* ix */) {
                 return detail::get_or_throw<int>("diag", "file_min_digits");
             },
             [](ImpactX & /* ix */, int const file_min_digits) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("file_min_digits", file_min_digits);
             },
             "The minimum number of digits (default: 6) used for the step\n"
             "number appended to the diagnostic file names."
         )
        .def_property("abort_on_warning_threshold",
             [](ImpactX & /* ix */){
                 return detail::get_or_throw<std::string>("impactx", "abort_on_warning_threshold");
             },
             [](ImpactX & ix, std::string const str_abort_on_warning_threshold) {
                 amrex::ParmParse pp_impactx("impactx");
                 pp_impactx.add("abort_on_warning_threshold", str_abort_on_warning_threshold);
                 // query input for warning logger variables and set up warning logger accordingly
                 ix.init_warning_logger();
             },
             "Set WarnPriority threshold to decide if ImpactX\n"
             "has to abort when a warning is recorded.\n"
             "Valid choices are: ['low', 'medium', 'high']."
        )
        .def_property("always_warn_immediately",
            [](ImpactX & /* ix */){
                 return detail::get_or_throw<int>("impactx", "always_warn_immediately");
              },
            [](ImpactX & /* ix */, int const always_warn_immediately) {
                amrex::ParmParse pp_impactx("impactx");
                pp_impactx.add("always_warn_immediately", always_warn_immediately);
            },
            "If set to 1, immediately prints every warning message\n"
            " as soon as it is generated."
        )
        // TODO this is an integer with 0 or 1 - can I just make this a boolean here?
        .def_property("abort_on_unused_inputs",
            [](ImpactX & /* ix */){
                return detail::get_or_throw<int>("amrex", "abort_on_unused_inputs");
            },
            [](ImpactX & ix, int const abort_on_unused_inputs) {
                amrex::ParmParse pp_amrex("amrex");
                pp_amrex.add("abort_on_unused_inputs", abort_on_unused_inputs);
                // query input for warning logger variables and set up warning logger accordingly
                ix.init_warning_logger();
            },
            "Configure simulation to abort AFTER it has run\n"
            "if there are unused parameters in the input."
        )

        .def("init_grids", &ImpactX::initGrids,
             "Initialize AMReX blocks/grids for domain decomposition & space charge mesh.\n\n"
             "This must come first, before particle beams and lattice elements are initialized."
        )
        .def("init_beam_distribution_from_inputs", &ImpactX::initBeamDistributionFromInputs)
        .def("init_lattice_elements_from_inputs", &ImpactX::initLatticeElementsFromInputs)
        .def("add_particles", &ImpactX::add_particles,
             py::arg("bunch_charge"),
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
        .def("reduced_beam_characteristics",
             [](ImpactX & ix) {
                return diagnostics::reduced_beam_characteristics(*ix.m_particle_container);
            },
            "Compute reduced beam characteristics like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters."
        )
        .def(
            "rho",
            [](ImpactX & ix, int const lev) { return &ix.m_rho.at(lev); },
            py::arg("lev"),
            py::return_value_policy::reference_internal,
            "charge density per level"
        )
        .def(
            "phi",
            [](ImpactX & ix, int const lev) { return &ix.m_phi.at(lev); },
            py::arg("lev"),
            py::return_value_policy::reference_internal,
            "scalar potential per level"
        )
        .def(
            "space_charge_field",
            [](ImpactX & ix, int const lev, std::string const comp) {
                return &ix.m_space_charge_field.at(lev).at(comp);
            },
            py::arg("lev"), py::arg("comp"),
            py::return_value_policy::reference_internal,
            "space charge force (vector: x,y,z) per level"
        )
        .def_readwrite("lattice",
            &ImpactX::m_lattice,
            "Access the accelerator element lattice."
        )
        .def_property("periods",
              [](ImpactX & /* ix */) {
                  return detail::get_or_throw<int>("lattice", "periods");
              },
              [](ImpactX & /* ix */, int periods) {
                  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(periods >= 1,
                                                   "lattice.periods must be >= 1");
                  amrex::ParmParse pp_lattice("lattice");
                  pp_lattice.add("periods", periods);
              },
              "The number of periods to repeat the lattice."
        )

        // from AmrCore->AmrMesh
        .def("Geom",
            //[](ImpactX const & ix, int const lev) { return ix.Geom(lev); },
            py::overload_cast< int >(&ImpactX::Geom, py::const_),
            py::arg("lev")
        )
        .def("DistributionMap",
            [](ImpactX const & ix, int const lev) { return ix.DistributionMap(lev); },
            //py::overload_cast< int >(&ImpactX::DistributionMap, py::const_),
            py::arg("lev")
        )
        .def("boxArray",
            [](ImpactX const & ix, int const lev) { return ix.boxArray(lev); },
            //py::overload_cast< int >(&ImpactX::boxArray, py::const_),
            py::arg("lev")
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
