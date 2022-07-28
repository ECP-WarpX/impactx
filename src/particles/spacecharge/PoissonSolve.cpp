/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "PoissonSolve.H"

#include <ablastr/fields/PoissonSolver.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_REAL.H>       // for ParticleReal

#include <cmath>


namespace impactx::spacecharge
{
    void PoissonSolve (
        ImpactXParticleContainer const & pc,
        std::unordered_map<int, amrex::MultiFab> & rho,
        std::unordered_map<int, amrex::MultiFab> & phi
    )
    {
        using namespace amrex::literals;

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
        amrex::ParticleReal const beta_s = std::sqrt(1.0_prt - 1.0_prt/std::pow(pt_ref, 2));
        // The beam particles and the corresponding box are all given in local coordinates
        // in which z is the direction of motion - this coincides with the direction of the momentum
        // of the reference particle.
        // After every T-to-Z transformation, Z aligns with the tangential vector of our reference
        // particle.
        std::array<amrex::Real, 3> const beta_xyz = {0.0, 0.0, beta_s};

        amrex::Real const mlmg_relative_tolerance = 1.e-7; // relative TODO: make smaller for SP
        amrex::Real const mlmg_absolute_tolerance = 0.0;   // ignored

        int const mlmg_max_iters = 100;
        int const mlmg_verbosity = 1;

        struct PoissonBoundaryHandler {
            amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const lobc = {
                amrex::LinOpBCType::Dirichlet,
                amrex::LinOpBCType::Dirichlet,
                amrex::LinOpBCType::Dirichlet
            };
            amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const hibc = {
                amrex::LinOpBCType::Dirichlet,
                amrex::LinOpBCType::Dirichlet,
                amrex::LinOpBCType::Dirichlet
            };
            //bool bcs_set = false;
            //std::array<bool, AMREX_SPACEDIM * 2> dirichlet_flag;
            //bool has_non_periodic = false;
        } poisson_boundary_handler;

        // create a vector to our fields, sorted by level
        amrex::Vector<amrex::MultiFab*> sorted_rho;
        amrex::Vector<amrex::MultiFab*> sorted_phi;

        for (int lev = 0; lev <= finest_level; ++lev) {
            sorted_rho.emplace_back(&rho[lev]);
            sorted_phi.emplace_back(&phi[lev]);
        }

        ablastr::fields::computePhi(
            sorted_rho,
            sorted_phi,
            beta_xyz,
            mlmg_relative_tolerance,
            mlmg_absolute_tolerance,
            mlmg_max_iters,
            mlmg_verbosity,
            pc.GetParGDB()->Geom(),
            pc.GetParGDB()->DistributionMap(),
            pc.GetParGDB()->boxArray(),
            poisson_boundary_handler
            /*
            do_single_precision_comms,
            this->ref_ratio,
            post_phi_calculation,
            gett_new(0),
            eb_farray_box_factory
            */
        );
    }
} // impactx::spacecharge
