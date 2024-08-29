/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <ImpactX.H>
#include <particles/transformation/CoordinateTransformation.H>

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
    auto get_or_throw (std::string const & prefix, std::string const & name)
    {
        T value;
        // TODO: if array do queryarr
        // bool const has_name = amrex::ParmParse(prefix).queryarr(name.c_str(), value);
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
            [](ImpactX const & /* ix */, std::string const & filename) {
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
                amrex::ParmParse const pp_amr("amr");
                pp_amr.getarr("n_cell", n_cell);
                return n_cell;
            },
            [](ImpactX & ix, std::array<int, AMREX_SPACEDIM> n_cell) {
                if (ix.initialized())
                    throw std::runtime_error("Read-only parameter after init_grids was called.");

                amrex::ParmParse pp_amr("amr");
                amrex::Vector<int> const n_cell_v(n_cell.begin(), n_cell.end());
                pp_amr.addarr("n_cell", n_cell_v);
            },
            "The number of grid points along each direction on the coarsest level."
        );

    for (int dir : {0, 1, 2}) {
        std::string const dir_str = std::vector<std::string>{"x", "y", "z"}.at(dir);
        std::string const bf_str = "blocking_factor_" + dir_str;
        impactx
            .def_property(bf_str.c_str(),
                  [bf_str](ImpactX & /* ix */) {
                      amrex::ParmParse pp_amr("amr");
                      std::vector<int> blocking_factor_dir;
                      pp_amr.queryarr(bf_str.c_str(), blocking_factor_dir);

                      return blocking_factor_dir;
                  },
                  [bf_str](ImpactX &ix, std::vector<int> blocking_factor_dir) {
                      if (ix.initialized())
                          throw std::runtime_error("Read-only parameter after init_grids was called.");

                      amrex::ParmParse pp_amr("amr");
                      pp_amr.addarr(bf_str.c_str(), blocking_factor_dir);
                  },
                  "AMReX blocking factor for a direction, per MR level."
            );
    }

    impactx
        // from amrex::AmrMesh
        .def_property("max_level",
            [](ImpactX & /* ix */){
                int max_level = 0;
                amrex::ParmParse pp_amr("amr");
                pp_amr.query("max_level", max_level);
                return max_level;
            },
            [](ImpactX & ix, int max_level) {
                if (ix.initialized())
                    throw std::runtime_error("Read-only parameter after init_grids was called.");

                amrex::ParmParse pp_amr("amr");
                pp_amr.add("max_level", max_level);
            },
            "The maximum mesh-refinement level for the simulation."
        )
        .def_property_readonly("finest_level",
            [](ImpactX & ix){ return ix.amr_data->finestLevel(); },
            "The currently finest level of mesh-refinement used. This is always less or equal to max_level."
        )

        .def_property("domain",
            [](ImpactX & /* ix */) {
                amrex::ParmParse const pp_geometry("geometry");
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
              [](ImpactX & /* ix */, std::vector<amrex::Real> frac) {
                  amrex::ParmParse pp_geometry("geometry");
                  pp_geometry.addarr("prob_relative", frac);
              },
              "The field mesh spans, per direction, multiple times the maximum physical extent of beam particles, as given by this factor."
        )

        .def_property("dynamic_size",
              [](ImpactX & /* ix */) {
                  amrex::ParmParse const pp_geometry("geometry");
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
            [](ImpactX & /* ix */, int const order) {
                amrex::ParmParse pp_ago("algo");
                pp_ago.add("particle_shape", order);
            },
            "Whether to calculate space charge effects."
        )
        .def_property("csr",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("algo", "csr");
            },
            [](ImpactX & /* ix */, bool const enable) {
                amrex::ParmParse pp_algo("algo");
                pp_algo.add("csr", enable);
            },
            "Enable or disable Coherent Synchrotron Radiation (CSR) calculations (default: disabled)."
        )
        .def_property("csr_bins",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("algo", "csr_bins");
            },
            [](ImpactX & /* ix */, int csr_bins) {
                amrex::ParmParse pp_algo("algo");
                pp_algo.add("csr_bins", csr_bins);
            },
            "Number of longitudinal bins used for CSR calculations (default: 150)."
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
        .def_property("poisson_solver",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<std::string>("algo", "poisson_solver");
            },
            [](ImpactX & /* ix */, std::string const poisson_solver) {
                if (poisson_solver != "multigrid" && poisson_solver != "fft") {
                    throw std::runtime_error("Poisson solver must be multigrid or fft but is: " + poisson_solver);
                }

                amrex::ParmParse pp_algo("algo");
                pp_algo.add("poisson_solver", poisson_solver);
            },
            "The numerical solver to solve the Poisson equation when calculating space charge effects. Either multigrid (default) or fft."
        )
        .def_property("mlmg_relative_tolerance",
              [](ImpactX & /* ix */) {
                  return detail::get_or_throw<bool>("algo", "mlmg_relative_tolerance");
              },
              [](ImpactX & /* ix */, amrex::Real const mlmg_relative_tolerance) {
                  amrex::ParmParse pp_algo("algo");
                  pp_algo.add("mlmg_relative_tolerance", mlmg_relative_tolerance);
              },
              "The relative precision with which the electrostatic space-charge fields should be calculated. "
              "More specifically, the space-charge fields are computed with an iterative Multi-Level Multi-Grid (MLMG) solver. "
              "This solver can fail to reach the default precision within a reasonable time."
        )
        .def_property("mlmg_absolute_tolerance",
              [](ImpactX & /* ix */) {
                  return detail::get_or_throw<bool>("algo", "mlmg_absolute_tolerance");
              },
              [](ImpactX & /* ix */, amrex::Real const mlmg_absolute_tolerance) {
                  amrex::ParmParse pp_algo("algo");
                  pp_algo.add("mlmg_absolute_tolerance", mlmg_absolute_tolerance);
              },
              "The absolute tolerance with which the space-charge fields should be calculated in units of V/m^2. "
              "More specifically, the acceptable residual with which the solution can be considered converged. "
              "In general this should be left as the default, but in cases where the simulation state changes very "
              "little between steps it can occur that the initial guess for the MLMG solver is so close to the "
              "converged value that it fails to improve that solution sufficiently to reach the "
              "mlmg_relative_tolerance value."
        )
        .def_property("mlmg_max_iters",
              [](ImpactX & /* ix */) {
                  return detail::get_or_throw<bool>("algo", "mlmg_max_iters");
              },
              [](ImpactX & /* ix */, int const mlmg_max_iters) {
                  amrex::ParmParse pp_algo("algo");
                  pp_algo.add("mlmg_max_iters", mlmg_max_iters);
              },
              "Maximum number of iterations used for MLMG solver for space-charge fields calculation. "
              "In case if MLMG converges but fails to reach the desired self_fields_required_precision, "
              "this parameter may be increased."
        )
        .def_property("mlmg_verbosity",
              [](ImpactX & /* ix */) {
                  return detail::get_or_throw<bool>("algo", "mlmg_verbosity");
              },
              [](ImpactX & /* ix */, int const mlmg_verbosity) {
                  amrex::ParmParse pp_algo("algo");
                  pp_algo.add("mlmg_verbosity", mlmg_verbosity);
              },
              "The verbosity used for MLMG solver for space-charge fields calculation. "
              "Currently MLMG solver looks for verbosity levels from 0-5. "
              "A higher number results in more verbose output."
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
        .def_property("particle_lost_diagnostics_backend",
                      [](ImpactX & /* ix */) {
                          return detail::get_or_throw<std::string>("diag", "backend");
                      },
                      [](ImpactX & /* ix */, std::string const backend) {
                          amrex::ParmParse pp_diag("diag");
                          pp_diag.add("backend", backend);
                      },
                      "Diagnostics for particles lost in apertures.\n\n"
                      "See the ``BeamMonitor`` element for backend values."
        )
        .def_property("abort_on_warning_threshold",
             [](ImpactX & /* ix */){
                 return detail::get_or_throw<std::string>("impactx", "abort_on_warning_threshold");
             },
             [](ImpactX & ix, std::string const & str_abort_on_warning_threshold) {
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

        .def_property("verbose",
            [](ImpactX & /* ix */){
                return detail::get_or_throw<int>("impactx", "verbose");
            },
            [](ImpactX & /* ix */, int const verbose) {
                amrex::ParmParse pp_impactx("impactx");
                pp_impactx.add("verbose", verbose);
            },
            "Controls how much information is printed to the terminal, when running ImpactX.\n"
            "``0`` for silent, higher is more verbose. Default is ``1``."
        )

        .def("deposit_charge",
            [](ImpactX & ix) {
                // transform from x',y',t to x,y,z
                transformation::CoordinateTransformation(
                        *(ix.amr_data)->m_particle_container,
                        CoordSystem::t);

                // Note: The following operation assume that
                // the particles are in x, y, z coordinates.

                // Resize the mesh, based on `m_particle_container` extent
                ix.ResizeMesh();

                // Redistribute particles in the new mesh in x, y, z
                ix.amr_data->m_particle_container->Redistribute();

                // charge deposition
                ix.amr_data->m_particle_container->DepositCharge(ix.amr_data->m_rho, ix.amr_data->refRatio());

                // transform from x,y,z to x',y',t
                transformation::CoordinateTransformation(*(ix.amr_data)->m_particle_container,
                                                         CoordSystem::s);
            },
            "Deposit charge in x,y,z."
        )

        .def("finalize", &ImpactX::finalize,
             "Deallocate all contexts and data."
        )
        .def("init_grids", &ImpactX::init_grids,
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
        // TODO: step
        .def("resize_mesh", &ImpactX::ResizeMesh,
             "Resize the mesh :py:attr:`~domain` based on the :py:attr:`~dynamic_size` and related parameters."
        )

        .def("particle_container",
             [](ImpactX & ix) -> ImpactXParticleContainer & {
                return *ix.amr_data->m_particle_container;
             },
             py::return_value_policy::reference_internal,
             "Access the beam particle container."
        )
        .def(
            "rho",
            [](ImpactX & ix, int const lev) { return &ix.amr_data->m_rho.at(lev); },
            py::arg("lev"),
            py::return_value_policy::reference_internal,
            "charge density per level"
        )
        .def(
            "phi",
            [](ImpactX & ix, int const lev) { return &ix.amr_data->m_phi.at(lev); },
            py::arg("lev"),
            py::return_value_policy::reference_internal,
            "scalar potential per level"
        )
        .def(
            "space_charge_field",
            [](ImpactX & ix, int lev, std::string const & comp) {
                return &ix.amr_data->m_space_charge_field.at(lev).at(comp);
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
            [](ImpactX const & ix, int const lev) { return ix.amr_data->Geom(lev); },
            py::arg("lev")
        )
        .def("DistributionMap",
            [](ImpactX const & ix, int const lev) { return ix.amr_data->DistributionMap(lev); },
            //py::overload_cast< int >(&ImpactX::DistributionMap, py::const_),
            py::arg("lev")
        )
        .def("boxArray",
            [](ImpactX const & ix, int const lev) { return ix.amr_data->boxArray(lev); },
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
            [](py::object const &){
#ifdef AMREX_USE_MPI
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "have_gpu",
            [](py::object const &){
#ifdef AMREX_USE_GPU
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "have_omp",
            [](py::object const &){
#ifdef AMREX_USE_OMP
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "gpu_backend",
            [](py::object const &){
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
