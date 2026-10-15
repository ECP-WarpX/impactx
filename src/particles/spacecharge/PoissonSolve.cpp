/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "PoissonSolve.H"

#include "initialization/Algorithms.H"
#include "particles/ChargeDeposition.H"

#include <ablastr/fields/MLMGOptions.H>
#include <ablastr/fields/PoissonSolver.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>       // for ParticleReal

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>


namespace impactx::particles::spacecharge
{
    void PoissonSolve (
        ImpactXParticleContainer const & pc,
        std::unordered_map<int, amrex::MultiFab> & rho,
        std::unordered_map<int, amrex::MultiFab> & rho_2d_out,
        std::unordered_map<int, amrex::MultiFab> & phi,
        amrex::Vector<amrex::IntVect> rel_ref_ratio,
        amrex::Vector<initialization::EBFactory const *> const & eb_factory
    )
    {
        BL_PROFILE("impactx::spacecharge::PoissonSolve");

        using namespace amrex::literals;
        using amrex::Math::powi;

        auto space_charge = get_space_charge_algo();

        // set space charge field to zero
        //   loop over refinement levels
        int const finest_level = phi.size() - 1u;
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::MultiFab &phi_at_level = phi.at(lev);
            // reset the values in phi to zero
            phi_at_level.setVal(0.);
        }

        // prepare parameters of the MLMG Poisson Solver
        //   relativistic beta=v/c of the reference particle
        amrex::ParticleReal const pt_ref = pc.GetRefParticle().pt;
        amrex::ParticleReal const beta_s = std::sqrt(1.0_prt - 1.0_prt/powi<2>(pt_ref));
        // The beam particles and the corresponding box are all given in local coordinates
        // in which z is the direction of motion - this coincides with the direction of the momentum
        // of the reference particle.
        // After every T-to-Z transformation, Z aligns with the tangential vector of our reference
        // particle.
        std::array<amrex::Real, 3> const beta_xyz = {0.0, 0.0, beta_s};

        amrex::ParmParse pp_algo("algo");
        std::string poisson_solver = "multigrid";
        pp_algo.queryAdd("poisson_solver", poisson_solver);
        const bool is_solver_igf_on_lev0 = poisson_solver == "fft";
        if (poisson_solver != "multigrid" && poisson_solver != "fft") {
            throw std::runtime_error("algo.poisson_solver must be multigrid or fft but is: " + poisson_solver);
        }
        if (space_charge == SpaceChargeAlgo::True_2D && poisson_solver != "fft") {
            throw std::runtime_error("algo.poisson_solver must be fft for SpaceChargeAlgo::True_2D");
        }
        if (space_charge == SpaceChargeAlgo::True_2p5D && poisson_solver != "fft") {
            throw std::runtime_error("algo.poisson_solver must be fft for SpaceChargeAlgo::True_2p5D");
        }
        bool const eb_enabled = !eb_factory.empty();
        if (eb_enabled &&
            (space_charge != SpaceChargeAlgo::True_3D || poisson_solver != "multigrid")) {
            throw std::runtime_error(
                "eb.shape requires algo.space_charge = 3D and algo.poisson_solver = multigrid");
        }

        // MLMG options
        //   Single precision: achievable relative residual is limited by float32
        //   round-off (machine epsilon ~1.2e-7).
#ifdef AMREX_USE_FLOAT
        amrex::Real mlmg_relative_tolerance = 1.e-4_rt; // relative (single precision)
#else
        amrex::Real mlmg_relative_tolerance = 1.e-7_rt; // relative (double precision)
#endif
        amrex::Real mlmg_absolute_tolerance = 0.0;   // ignored
        pp_algo.queryAddWithParser("mlmg_relative_tolerance", mlmg_relative_tolerance);
        pp_algo.queryAddWithParser("mlmg_absolute_tolerance", mlmg_absolute_tolerance);

        int mlmg_max_iters = 100;
        int mlmg_verbosity = 1;
        pp_algo.queryAddWithParser("mlmg_max_iters", mlmg_max_iters);
        pp_algo.queryAddWithParser("mlmg_verbosity", mlmg_verbosity);

        ablastr::fields::MLMGOptions mlmg_options;
        mlmg_options.relative_tolerance = mlmg_relative_tolerance;
        mlmg_options.absolute_tolerance = mlmg_absolute_tolerance;
        mlmg_options.max_iters = mlmg_max_iters;
        mlmg_options.verbosity = mlmg_verbosity;

        // flatten rho to 2D; store it in the output so it can be accessed after
        // the solve (e.g. via sim.rho), like the solved potential phi
        if (space_charge == SpaceChargeAlgo::True_2D || space_charge == SpaceChargeAlgo::True_2p5D) {
            auto geom_3d = pc.GetParGDB()->Geom();
            amrex::Box domain_3d = geom_3d[0].Domain();  // whole simulation index space (level 0)
            rho_2d_out = project_charge_to_2D(rho, domain_3d);
        }

        // create a vector to our fields, sorted by level
        amrex::Vector<amrex::MultiFab*> sorted_rho;
        amrex::Vector<amrex::MultiFab*> sorted_phi;

        amrex::Vector<amrex::MultiFab> phi_2d(finest_level+1);

        // create phi_2d and sort rho/phi pointers
        for (int lev = 0; lev <= finest_level; ++lev) {
            if (space_charge == SpaceChargeAlgo::True_2D || space_charge == SpaceChargeAlgo::True_2p5D) {
                int nz = pc.GetParGDB()->Geom(lev).Domain().length(2);
                if (nz == 1) {
                    sorted_phi.emplace_back(&phi[lev]);
                } else {
                    // 2D phi
                    auto & r2d = rho_2d_out[lev];
                    auto nGrow = phi[lev].nGrowVect();
                    nGrow[2] = 0;
                    phi_2d[lev].define(r2d.boxArray(), r2d.DistributionMap(), r2d.nComp(), nGrow);
                    sorted_phi.emplace_back(&phi_2d[lev]);
                }

                sorted_rho.emplace_back(&rho_2d_out[lev]);
            }
            else if (space_charge == SpaceChargeAlgo::True_3D) {
                sorted_rho.emplace_back(&rho[lev]);
                sorted_phi.emplace_back(&phi[lev]);
            }
        }

        const bool is_igf_2d = (space_charge == SpaceChargeAlgo::True_2D || space_charge == SpaceChargeAlgo::True_2p5D);
        const bool do_single_precision_comms = false;
        std::optional<amrex::Vector<initialization::EBFactory const *>> eb_farray_box_factory;
        if (eb_enabled) {
            eb_farray_box_factory = eb_factory;
        }
        // EB: no boundary handler, i.e., phi = 0 on the EB and on all domain faces
        ablastr::fields::computePhi(
            sorted_rho,
            sorted_phi,
            beta_xyz,
            mlmg_options,
            pc.GetParGDB()->Geom(),
            pc.GetParGDB()->DistributionMap(),
            pc.GetParGDB()->boxArray(),
            ablastr::utils::enums::GridType::Collocated,
            is_solver_igf_on_lev0,
            is_igf_2d,
            eb_enabled,
            do_single_precision_comms,
            rel_ref_ratio,
            std::nullopt,  // post_phi_calculation
            std::nullopt,  // boundary_handler
            std::nullopt,  // current_time
            eb_farray_box_factory
        );

        // We may need to copy phi from phi_2d
        if (space_charge == SpaceChargeAlgo::True_2D || space_charge == SpaceChargeAlgo::True_2p5D) {
            for (int lev=0; lev<=finest_level; lev++) {
                if (&(phi[lev]) != sorted_phi[lev]) {
                    phi[lev].ParallelCopy(*sorted_phi[lev]);
                }
            }
        }

        // fill boundary
        for (int lev=0; lev<=finest_level; lev++)
        {
            amrex::MultiFab & phi_at_level = phi.at(lev);
            phi_at_level.FillBoundary(pc.GetParGDB()->Geom()[lev].periodicity());
        }
    }
}  // impactx::particles::spacecharge
