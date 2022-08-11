/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "ForceFromSelfFields.H"

#include <ablastr/fields/PoissonSolver.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>       // for Real
#include <AMReX_SPACE.H>      // for AMREX_D_DECL


namespace impactx::spacecharge
{
    void ForceFromSelfFields (
        std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > & space_charge_force,
        std::unordered_map<int, amrex::MultiFab> const & phi,
        amrex::Vector<amrex::Geometry> const & geom
    )
    {
        BL_PROFILE("impactx::spacecharge::ForceFromSelfFields");

        using namespace amrex::literals;

        // loop over refinement levels
        int const finest_level = phi.size() - 1u;
        for (int lev = 0; lev <= finest_level; ++lev) {

            // stencil coefficients: 0.5 * inverse cell size
            auto const &gm = geom[lev];
            auto const dr = gm.CellSizeArray();
            amrex::GpuArray<amrex::Real, 3> const inv2dr{AMREX_D_DECL(0.5_rt/dr[0], 0.5_rt/dr[1], 0.5_rt/dr[2])};

            // reset the values in space_charge_force to zero
            space_charge_force.at(lev).at("x").setVal(0.);
            space_charge_force.at(lev).at("y").setVal(0.);
            space_charge_force.at(lev).at("z").setVal(0.);

            for (amrex::MFIter mfi(phi.at(lev)); mfi.isValid(); ++mfi) {

                amrex::Box bx = mfi.validbox();
                auto const phi_arr = (phi.at(lev))[mfi].const_array();

                auto scf_arr_x = space_charge_force[lev]["x"][mfi].array();
                auto scf_arr_y = space_charge_force[lev]["y"][mfi].array();
                auto scf_arr_z = space_charge_force[lev]["z"][mfi].array();

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    scf_arr_x(i, j, k) = inv2dr[0] * (phi_arr(i-1, j, k) - phi_arr(i+1, j, k));
                    scf_arr_y(i, j, k) = inv2dr[1] * (phi_arr(i, j-1, k) - phi_arr(i, j+1, k));
                    scf_arr_z(i, j, k) = inv2dr[2] * (phi_arr(i, j, k-1) - phi_arr(i, j, k+1));
                });
            }
        }
    }
} // namespace impactx::spacecharge
