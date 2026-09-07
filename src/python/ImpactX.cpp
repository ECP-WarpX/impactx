/* Copyright 2022-2026 The ImpactX Community
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <ImpactX.H>
#include <diagnostics/FilePrefix.H>
#include <initialization/Algorithms.H>
#include <initialization/InitDistribution.H>
#include <particles/ChargeDeposition.H>
#include <particles/CovarianceMatrix.H>
#include <particles/transformation/CoordinateTransformation.H>
#include <particles/ParticleBoundary.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#if defined(AMREX_DEBUG) || defined(DEBUG)
#   include <cstdio>
#endif
#include <optional>
#include <string>
#include <type_traits>
#include <variant>


namespace py = pybind11;
using namespace impactx;

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
        using V = std::decay_t<T>;
        V value;

        bool has_name = false;
        // TODO: if array do queryarr
        // has_name = amrex::ParmParse(prefix).queryarr(name.c_str(), value);
        if constexpr (std::is_same_v<V, bool> || std::is_same_v<V, std::string>)
            has_name = amrex::ParmParse(prefix).query(name.c_str(), value);
        else
            has_name = amrex::ParmParse(prefix).queryWithParser(name.c_str(), value);

        if (!has_name)
            throw std::runtime_error(prefix + "." + name + " is not set yet");
        return value;
    }
}

void init_ImpactX (py::module& m)
{
    py::class_<ImpactX> impactx(m, "ImpactX", py::dynamic_attr());

    using ImpactXHooks = std::unordered_map<std::string, std::function<void(ImpactX *)> >;
    py::bind_map<ImpactXHooks>(m, "UnorderedMap");

    impactx
        .def(py::init<>())

        .def("load_inputs_file",
            [](ImpactX const & /* ix */, std::string const & filename) {
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
        std::string const mgs_str = "max_grid_size_" + dir_str;
        impactx
            .def_property(bf_str.c_str(),
                  [bf_str](ImpactX & /* ix */) {
                      amrex::ParmParse pp_amr("amr");
                      std::vector<int> blocking_factor_dir;
                      pp_amr.queryarr(bf_str.c_str(), blocking_factor_dir);

                      return blocking_factor_dir;
                  },
                  [bf_str, mgs_str](ImpactX &ix, std::vector<int> blocking_factor_dir) {
                      if (ix.initialized())
                          throw std::runtime_error("Read-only parameter after init_grids was called.");

                      amrex::ParmParse pp_amr("amr");
                      pp_amr.addarr(bf_str.c_str(), blocking_factor_dir);
                      pp_amr.addarr(mgs_str.c_str(), blocking_factor_dir);
                  },
                  "AMReX blocking factor for a direction, per MR level. "
                  "Also sets max_grid_size_x/_y/_z in this direction to the "
                  "same value; assign that property afterwards to overwrite. "
                  "The all-direction max_grid_size does not overwrite it, "
                  "since per-direction values always take precedence."
            )
            .def_property(mgs_str.c_str(),
                  [mgs_str](ImpactX & /* ix */) {
                      amrex::ParmParse pp_amr("amr");
                      std::vector<int> max_grid_size_dir;
                      pp_amr.queryarr(mgs_str.c_str(), max_grid_size_dir);

                      return max_grid_size_dir;
                  },
                  [mgs_str](ImpactX &ix, std::vector<int> max_grid_size_dir) {
                      if (ix.initialized())
                          throw std::runtime_error("Read-only parameter after init_grids was called.");

                      amrex::ParmParse pp_amr("amr");
                      pp_amr.addarr(mgs_str.c_str(), max_grid_size_dir);
                  },
                  "AMReX maximum box size for a direction, per MR level. "
                  "Boxes larger than this are chopped, aim for one box per "
                  "parallel process. Takes precedence over the all-direction "
                  "max_grid_size, no matter in which order they are assigned."
            );
    }

    impactx
        .def_property("max_grid_size",
            [](ImpactX & /* ix */) {
                amrex::ParmParse pp_amr("amr");
                std::vector<int> max_grid_size;
                pp_amr.queryarr("max_grid_size", max_grid_size);

                return max_grid_size;
            },
            [](ImpactX &ix, std::vector<int> max_grid_size) {
                if (ix.initialized())
                    throw std::runtime_error("Read-only parameter after init_grids was called.");

                amrex::ParmParse pp_amr("amr");
                pp_amr.addarr("max_grid_size", max_grid_size);
            },
            "AMReX maximum box size along all directions, per MR level. "
            "Boxes larger than this are chopped; aim for one box per "
            "parallel process. Per-direction values (max_grid_size_x/_y/_z, "
            "also set by blocking_factor_x/_y/_z) take precedence, no matter "
            "in which order they are assigned; assign those to change a "
            "direction that is already set."
        );

    impactx
        // from amrex::AmrMesh
        .def_property("max_level",
            [](ImpactX & /* ix */){
                int max_level = 0;
                amrex::ParmParse pp_amr("amr");
                pp_amr.queryWithParser("max_level", max_level);
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

                if (ix.initialized())
                    ix.ResizeMesh();
            },
            "The physical extent of the full simulation domain, relative to the reference particle position, in meters."
        )

        .def_property("prob_relative",
              [](ImpactX & /* ix */) {
                  amrex::ParmParse const pp_geometry("geometry");
                  std::vector<amrex::Real> frac;
                  pp_geometry.getarr("prob_relative", frac);
                  return frac;
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
        .def_property("strang_split",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("algo", "strang_split");
            },
            [](ImpactX & /* ix */, bool const enable) {
                amrex::ParmParse pp_algo("algo");
                pp_algo.add("strang_split", enable);
            },
            "Compose the collective effect kicks and the element transport in a "
            "second-order, time-symmetric Strang split (default: enabled) or to first order."
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
        .def_property("isr",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("algo", "isr");
            },
            [](ImpactX & /* ix */, bool const enable) {
                amrex::ParmParse pp_algo("algo");
                pp_algo.add("isr", enable);
            },
            "Enable or disable Incoherent Synchrotron Radiation (ISR) calculations (default: disabled)."
        )
        .def_property("isr_order",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("algo", "isr_order");
            },
            [](ImpactX & /* ix */, int isr_order) {
                amrex::ParmParse pp_algo("algo");
                pp_algo.add("isr_order", isr_order);
            },
            "Number of terms in the Taylor series retained for quantum effects (default: 1)."
        )
        .def_property("isr_on_ref_part",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("algo", "isr_on_ref_part");
            },
            [](ImpactX & /* ix */, bool isr_on_ref_part) {
                amrex::ParmParse pp_algo("algo");
                pp_algo.add("isr_on_ref_part", isr_on_ref_part);
            },
            "Flag to determine whether ISR radiation loss is applied to the reference particle (default: False)."
        )
        .def_property("spin",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("algo", "spin");
            },
            [](ImpactX & /* ix */, bool const enable) {
                amrex::ParmParse pp_algo("algo");
                pp_algo.add("spin", enable);
            },
            "Enable or disable particle spin tracking (default: disabled)."
        )
        .def_property("eigenemittances",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<bool>("diag", "eigenemittances");
            },
            [](ImpactX & /* ix */, bool const enable) {
                amrex::ParmParse pp_diag("diag");
                pp_diag.add("eigenemittances", enable);
            },
            "Enable or disable eigenemittance diagnostic calculations (default: disabled)."
        )
        .def_property("space_charge",
            [](ImpactX & /* ix */) -> std::string {
                return detail::get_or_throw<std::string>("algo", "space_charge");
            },
            [](ImpactX & /* ix */, std::variant<bool, std::string> space_charge_v) {
                if (std::holds_alternative<bool>(space_charge_v)) {
                    amrex::ParmParse pp_algo("algo");
                    if (std::get<bool>(space_charge_v)) {
                        // TODO: boolean True is deprecated since 25.03, remove some time after
                        py::warnings::warn(
                            "sim.space_charge = True is deprecated, please use space_charge = \"3D\"",
                            PyExc_DeprecationWarning,
                            2
                        );
                        pp_algo.add("space_charge", std::string("3D"));
                    } else {
                        // map boolean False to "false" / off
                        pp_algo.add("space_charge", std::string("false"));
                    }
                }
                else
                {
                    std::string const space_charge = std::get<std::string>(space_charge_v);
                    if (space_charge != "false" && space_charge != "off" && space_charge != "2D" && space_charge != "3D" && space_charge != "Gauss3D" && space_charge != "Gauss2p5D"  && space_charge != "2p5D") {
                        throw std::runtime_error("Space charge model must be 2D, 3D, Gauss3D, Gauss2p5D, or 2p5D but is: " + space_charge);
                    }
                    amrex::ParmParse pp_algo("algo");
                    pp_algo.add("space_charge", space_charge);
                }
            },
            "The model to be used when calculating space charge effects. Either off, 2D, 3D, Gauss3D, Gauss2p5D, or 2p5D."
        )
        .def_property("space_charge_gauss_nint",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<int>("algo.space_charge", "gauss_nint");
            },
            [](ImpactX & /* ix */, int const gauss_nint) {
                if (gauss_nint < 1) {
                    throw std::runtime_error("space_charge_gauss_nint must be strictly positive");
                }

                amrex::ParmParse pp_algo("algo.space_charge");
                pp_algo.add("gauss_nint", gauss_nint);
            },
            "Number of steps for computing the integrals (default: ``101``)."
        )
        .def_property("space_charge_gauss_charge_z_bins",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<int>("algo.space_charge", "gauss_charge_z_bins");
            },
            [](ImpactX & /* ix */, int const gauss_charge_z_bins) {
                if (gauss_charge_z_bins < 1) {
                    throw std::runtime_error("space_charge_gauss_charge_z_bins must be strictly positive");
                }

                amrex::ParmParse pp_algo("algo.space_charge");
                pp_algo.add("gauss_charge_z_bins", gauss_charge_z_bins);
            },
            "Number of longitudinal bins for computing the linear charge density (default: ``129``)."
        )
        .def_property("space_charge_gauss_taylor_delta",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<int>("algo.space_charge", "gauss_taylor_delta");
            },
            [](ImpactX & /* ix */, amrex::Real const gauss_taylor_delta) {
                if (gauss_taylor_delta < 0) {
                    throw std::runtime_error("space_charge_gauss_taylor_delta must be strictly positive");
                }
                if (gauss_taylor_delta > 0.05) {
                    throw std::runtime_error("space_charge_gauss_taylor_delta must be less than 0.05");
                }
                amrex::ParmParse pp_algo("algo.space_charge");
                pp_algo.add("gauss_taylor_delta", gauss_taylor_delta);
            },
            "Initial region for computing the integrals (default: ``0.01``)."
        )
        .def_property("space_charge_gauss_long_scale",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<amrex::Real>("algo.space_charge", "gauss_long_scale");
            },
            [](ImpactX & /* ix */, amrex::Real const gauss_long_scale) {
                if (gauss_long_scale < 0) {
                    throw std::runtime_error("space_charge_gauss_long_scale must be strictly positive");
                }
                amrex::ParmParse pp_algo("algo.space_charge");
                pp_algo.add("gauss_long_scale", gauss_long_scale);
            },
            "Longitudinal space charge scale for the Gauss2p5D space charge model. "
            "Approximation affecting only the longitudinal momentum (``pt``) kick. "
            "If not set, it defaults to ``1.103 * gamma * sigma_z``, estimated in-situ from the current "
            "reduced beam characteristics, which is a typical value when comparing to a 3D model."
        )
        .def_property("space_charge_num_longitudinal_bins",
            [](ImpactX & /* ix */) {
                return detail::get_or_throw<int>("algo.space_charge", "num_longitudinal_bins");
            },
            [](ImpactX & /* ix */, int const num_longitudinal_bins) {
                if (num_longitudinal_bins < 1) {
                    throw std::runtime_error("space_charge_num_longitudinal_bins must be strictly positive");
                }
                amrex::ParmParse pp_algo("algo.space_charge");
                pp_algo.add("num_longitudinal_bins", num_longitudinal_bins);
            },
            "Number of longitudinal bins for 2.5D space charge calculation (default: ``100``)."
        )
        .def_property("space_charge_apply_longitudinal_kick",
             [](ImpactX & /* ix */) {
                 return detail::get_or_throw<bool>("algo.space_charge", "apply_longitudinal_kick");
             },
             [](ImpactX & /* ix */, bool const enable) {
                 amrex::ParmParse pp_algo("algo.space_charge");
                 pp_algo.add("apply_longitudinal_kick", enable);
             },
             "Enable or disable longitudinal space charge kick in 2.5D space charge solver (default: enabled).\n"
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
        .def_property("particle_bc",
            [](ImpactX & /* ix */) -> std::string {
                return amrex::getEnumNameString(particles::get_particle_boundary_condition());
            },
            [](ImpactX & /* ix */, std::string const particle_bc) {
                auto const valid_names = amrex::getEnumNameStrings<particles::ParticleBC>();
                if (std::find(valid_names.begin(), valid_names.end(), particle_bc) == valid_names.end()) {
                    std::string msg = "Particle boundary condition must be one of: ";
                    for (auto const& name : valid_names) {
                        msg += name + ", ";
                    }
                    msg.erase(msg.size() - 2);
                    msg += " but is: " + particle_bc;
                    throw std::runtime_error(msg);
                }

                amrex::ParmParse("algo").add("particle_bc", particle_bc);
            },
            "Optional methods to apply a longitudinal particle boundary condition."
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
        .def_property("diag_file_prefix",
             [](ImpactX & /* ix */) {
                 return diagnostics::FilePrefix();
             },
             [](ImpactX & /* ix */, std::string const & file_prefix) {
                 amrex::ParmParse pp_diag("diag");
                 pp_diag.add("file_prefix", file_prefix);
             },
             "Root directory for diagnostic output (default: ``diags``)."
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
                 //   the warning logger needs AMReX: if it is not yet initialized,
                 //   ImpactX::init_grids sets the logger up from the same input
                 if (amrex::Initialized()) {
                     ix.init_warning_logger();
                 }
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
                //   the warning logger needs AMReX: if it is not yet initialized,
                //   ImpactX::init_grids sets the logger up from the same input
                if (amrex::Initialized()) {
                    ix.init_warning_logger();
                }
            },
            "Configure simulation to abort AFTER it has run\n"
            "if there are unused parameters in the input."
        )

        .def_property("omp_threads",
            [](ImpactX & /* ix */){
                return detail::get_or_throw<std::string>("amrex", "omp_threads");
            },
            [](ImpactX & /* ix */, std::variant<int, std::string> omp_threads_var) {
                std::visit([&]( auto && omp_threads) {
                    amrex::ParmParse pp_amrex("amrex");
                pp_amrex.add("omp_threads", omp_threads);
                }, omp_threads_var);
            },
            "Controls the number of OpenMP threads to use (ImpactX default: \"nosmt\").\n"
            "https://amrex-codes.github.io/amrex/docs_html/InputsComputeBackends.html."
        )

        .def_property("verbose",
            [](ImpactX & /* ix */){
                return detail::get_or_throw<int>("impactx", "verbose");
            },
            [](ImpactX & /* ix */, int const verbose) {
                amrex::ParmParse pp_impactx("impactx");
                pp_impactx.add("verbose", verbose);
                // AMReX init/finalize
                amrex::ParmParse pp_amrex("amrex");
                pp_amrex.add("verbose", verbose);
            },
            "Controls how much information is printed to the terminal, when running ImpactX.\n"
            "``0`` for silent, higher is more verbose. Default is ``1``."
        )

        .def_property("tiny_profiler",
            [](ImpactX & /* ix */){
                return detail::get_or_throw<bool>("tiny_profiler", "enabled");
            },
            [](ImpactX & /* ix */, bool enabled) {
                amrex::ParmParse pp_tiny_profiler("tiny_profiler");
                pp_tiny_profiler.add("enabled", enabled);
            },
            "This parameter can be used to disable tiny profiling including CArena memory profiling at runtime."
        )
        .def_property("memory_profiler",
            [](ImpactX & /* ix */){
                return detail::get_or_throw<bool>("tiny_profiler", "memprof_enabled");
            },
            [](ImpactX & /* ix */, bool memprof_enabled) {
                amrex::ParmParse pp_tiny_profiler("tiny_profiler");
                pp_tiny_profiler.add("memprof_enabled", memprof_enabled);
            },
            "This parameter can be used to disable tiny profiler's memory arena profiling at runtime. "
            "If tiny_profiler.enabled is false, this parameter has no effects."
        )
        .def_property("tiny_profiler_file",
            [](ImpactX & /* ix */){
                return detail::get_or_throw<std::string>("tiny_profiler", "output_file");
            },
            [](ImpactX & /* ix */, std::string output_file) {
                amrex::ParmParse pp_tiny_profiler("tiny_profiler");
                pp_tiny_profiler.add("output_file", output_file);
            },
            "If this parameter is empty, the output of tiny profiling is dumped on the default out stream of AMReX. "
            "If it's not empty, it specifies the file name for the output. "
            "Note that /dev/null is a special name that mean a null file."
        )
        .def_readwrite("hook",
            &ImpactX::m_hook,
            "User-defined function hooks that are called, e.g, during tracking."
        )

        .def("deposit_charge",
            [](ImpactX & ix) {
                // transform from x',y',t to x,y,z
                transformation::CoordinateTransformation(
                        *(ix.amr_data)->track_particles.m_particle_container,
                        CoordSystem::t);

                // Note: The following operation assume that
                // the particles are in x, y, z coordinates.

                // Resize the mesh, based on `m_particle_container` extent
                ix.ResizeMesh();

                // Redistribute particles in the new mesh in x, y, z
                ix.amr_data->track_particles.m_particle_container->Redistribute();

                // charge deposition
                ix.amr_data->track_particles.m_particle_container->DepositCharge(
                    ix.amr_data->track_particles.m_rho,
                    ix.amr_data->refRatio());

                // for the 2D models, also produce the x-y projected (charge-
                // conserving) density, so it is available via sim.rho without a
                // Poisson solve
                auto const space_charge = get_space_charge_algo();
                if (space_charge == SpaceChargeAlgo::True_2D ||
                    space_charge == SpaceChargeAlgo::True_2p5D)
                {
                    auto geom_3d = ix.amr_data->track_particles
                        .m_particle_container->GetParGDB()->Geom();
                    amrex::Box const domain_3d = geom_3d[0].Domain();
                    ix.amr_data->track_particles.m_rho_2d = project_charge_to_2D(
                        ix.amr_data->track_particles.m_rho, domain_3d);
                }

                // transform from x,y,z to x',y',t
                transformation::CoordinateTransformation(
                    *(ix.amr_data)->track_particles.m_particle_container,
            CoordSystem::s
                );
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
        .def("init_envelope",
            [](ImpactX & ix, RefPart ref, distribution::KnownDistributions distr, std::optional<amrex::Real> intensity) {
                ix.amr_data->track_envelope.m_ref = ref;
                ix.amr_data->track_envelope.m_env = initialization::create_envelope(distr, intensity);
            },
            py::arg("ref"), py::arg("distr"), py::arg("intensity") = py::none(),
            "Envelope tracking mode:"
            "Create a 6x6 covariance matrix from a distribution and then initialize "
            "the simulation for envelope tracking relative to a reference particle."
        )
        .def("add_particles", &ImpactX::add_particles,
             py::arg("bunch_charge"),
             py::arg("distr"), py::arg("npart"),
             py::arg("spin_distr") = py::none(),
             "Particle tracking mode:"
             "Generate and add n particles to the particle container.\n\n"
             "Will also resize the geometry based on the updated particle\n"
             "distribution's extent and then redistribute particles in according\n"
             "AMReX grid boxes."
        )

        .def("evolve",  /** TODO: deprecated API. Only for internal use. Remove after a few releases. */
             [](ImpactX & ix) {
                py::warnings::warn(
                    "Warning: evolve() is deprecated and will soon be removed. "
                    "Use track_particles() instead.",
                    PyExc_DeprecationWarning,
                    2
                );
                ix.evolve();
             },
             "Run the main simulation loop."
        )
        .def("track_particles", &ImpactX::track_particles,
             "Run the particle tracking simulation loop."
        )
        .def("track_envelope", &ImpactX::track_envelope,
             "Run the envelope tracking simulation loop."
        )
        .def("track_reference", &ImpactX::track_reference,
             "Run the reference orbit tracking simulation loop."
        )


        .def_property_readonly("tracking_period",
            [](ImpactX & ix) { return ix.m_tracking_state.m_period; },
            "For tracking hooks/callbacks, the period in the lattice (e.g., turn or channel period)"
        )
        .def_property_readonly("tracking_step",
            [](ImpactX & ix) { return ix.m_tracking_state.m_step; },
            "For tracking hooks/callbacks, a global step of the simulation.\n\n"
            "A state of internal simulation steps, increments also for space charge slice steps in elements.\n"
            "We start in 'step 0' (initial state)."
        )
        .def_property_readonly("tracking_element",
            [](ImpactX & ix) { return ix.m_tracking_state.m_element; },
            "For tracking hooks/callbacks, the current lattice element."
        )

        .def("resize_mesh", &ImpactX::ResizeMesh,
             "Resize the mesh :py:attr:`~domain` based on the :py:attr:`~dynamic_size` and related parameters."
        )

        .def("particle_container",
             [](ImpactX & ix) -> ImpactXParticleContainer & {
                py::warnings::warn(
                    "particle_container() is deprecated. Use sim.beam instead.",
                    PyExc_DeprecationWarning,
                    2
                );
                return *ix.amr_data->track_particles.m_particle_container;
             },
             py::return_value_policy::reference_internal,
             "Access the beam particle container.\n\n"
             "Deprecated: use ``sim.beam``."
        )
        // Getter-only property is intentional: it returns the live mutable
        // ImpactXParticleContainer by reference.
        // A writable property (with assignment) would be ambiguous:
        // alias (Pythonic) or copy-in (safe) for the simulation-owned particle container.
        .def_property_readonly("beam",
            [](ImpactX & ix) -> ImpactXParticleContainer & {
                return *ix.amr_data->track_particles.m_particle_container;
            },
            "Access the beam particle container."
        )
        .def_property_readonly("envelope",
            [](ImpactX & ix) -> Envelope & {
                if (!ix.amr_data->track_envelope.m_env.has_value()) {
                    throw std::runtime_error(
                        "sim.envelope is only available after init_envelope()."
                    );
                }
                return ix.amr_data->track_envelope.m_env.value();
            },
            "Access the beam envelope (6x6 covariance matrix and beam intensity)\n"
            "used for envelope tracking. Only available after init_envelope()."
        )
        .def(
            "rho",
            [](py::object self, int const lev) -> py::object {
                ImpactX & ix = self.cast<ImpactX &>();
                auto & tp = ix.amr_data->track_particles;
                auto const space_charge = get_space_charge_algo();
                // The 2D models solve the x-y projected charge density: return
                // that (stored like the solved potential phi, populated by the
                // last deposit / Poisson solve); the 3D charge otherwise.
                auto & rho_per_level =
                    (space_charge == SpaceChargeAlgo::True_2D ||
                     space_charge == SpaceChargeAlgo::True_2p5D)
                    ? tp.m_rho_2d : tp.m_rho;
                return py::cast(
                    &rho_per_level.at(lev),
                    py::return_value_policy::reference_internal, self);
            },
            py::arg("lev"),
            "charge density per level.\n\n"
            "For the 2D space-charge models (``2D`` / ``2.5D``) this returns the\n"
            "x-y projected (2D) charge density that is actually solved; for the\n"
            "3D models it returns the 3D charge density."
        )
        .def(
            "phi",
            [](ImpactX & ix, int const lev) { return &ix.amr_data->track_particles.m_phi.at(lev); },
            py::arg("lev"),
            py::return_value_policy::reference_internal,
            "scalar potential per level"
        )
        .def(
            "space_charge_field",
            [](ImpactX & ix, int lev, std::string const & comp) {
                return &ix.amr_data->track_particles.m_space_charge_field.at(lev).at(comp);
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
}
