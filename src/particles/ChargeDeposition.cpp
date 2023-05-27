/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"

#include <ablastr/coarsen/average.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/particles/DepositCharge.H>

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParticleTile.H>

#include <array>


namespace impactx
{
    void
    ImpactXParticleContainer::DepositCharge (
        std::unordered_map<int, amrex::MultiFab> & rho,
        amrex::Vector<amrex::IntVect> const & ref_ratio)
    {
        // reset the values in rho to zero
        int const nLevel = this->finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            rho.at(lev).setVal(0.);
        }

        // loop fine-to-coarse over refinement levels
        for (int lev = nLevel; lev >= 0; --lev) {
            amrex::MultiFab & rho_at_level = rho.at(lev);

            // get simulation geometry information
            amrex::Geometry const & gm = this->Geom(lev);

            // Loop over particle tiles and deposit charge on each level
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                amrex::FArrayBox local_rho_fab;

                using ParIt = ImpactXParticleContainer::iterator;
                for (ParIt pti(*this, lev); pti.isValid(); ++pti) {
                    // preparing access to particle data: SoA of Reals
                    auto & AMREX_RESTRICT soa_real = pti.GetStructOfArrays().GetRealData();
                    // after https://github.com/ECP-WarpX/WarpX/pull/2838 add const:
                    auto const wp = soa_real[RealSoA::w];
                    int const * const AMREX_RESTRICT ion_lev = nullptr;

                    // physical lower corner of the current box
                    //   Note that this includes guard cells since it is after tilebox.grow
                    amrex::Box tilebox = pti.tilebox();
                    tilebox.grow(rho_at_level.nGrowVect());
                    amrex::RealBox const grid_box{tilebox, gm.CellSize(), gm.ProbLo()};
                    amrex::Real const * const AMREX_RESTRICT xyzmin_ptr = grid_box.lo();
                    std::array<amrex::Real, 3> const xyzmin = {xyzmin_ptr[0], xyzmin_ptr[1], xyzmin_ptr[2]};

                    // mesh-refinement: for when we do not deposit on the same level
                    // note: would need to communicate the deposited-to boxes afterwards
                    //int const depos_lev = lev;
                    // mesh refinement ratio between lev and depos_lev
                    //auto const rel_ref_ratio = ref_ratio.at(depos_lev) / ref_ratio.at(lev);
                    amrex::ignore_unused(ref_ratio);

                    // in SI [C]
                    amrex::ParticleReal const charge = m_refpart.charge;

                    // cell size of the mesh to deposit to
                    std::array<amrex::Real, 3> const & AMREX_RESTRICT dx = {gm.CellSize(0), gm.CellSize(1), gm.CellSize(2)};

                    // RZ modes (unused)
                    int const n_rz_azimuthal_modes = 0;

                    ablastr::particles::deposit_charge<ImpactXParticleContainer>
                            (pti, wp, charge, ion_lev, &rho_at_level,
                             local_rho_fab,
                             m_particle_shape.value(),
                             dx, xyzmin, n_rz_azimuthal_modes);
                }
            }

            // TODO: Call portion's of WarpX' SyncRho from fine to coarser levels
            //   TODO: do coarsening to a temp, local lev-1 patch
            //   TODO: start communicating this patch into rho_at_level_minus_1
            //   note: this can either move parts from WarpXComm.cpp (SyncRho & AddRhoFromFineLevelandSumBoundary)
            //         or use code from SyncPhi in ABLASTR's PoissonSolve (in opposite order: FP->CP here)
            // needed for solving the levels by levels:
            // - coarser level is initial guess for finer level
            // - coarser level provides boundary values for finer level patch
            // Interpolation from phi[lev] to phi[lev+1]
            // (This provides both the boundary conditions and initial guess for phi[lev+1])

            if (lev > 0) {
                // Allocate rho_cp with the same distribution map as lev
                amrex::BoxArray ba = rho[lev].boxArray();
                const amrex::IntVect &refratio = ref_ratio[lev - 1];
                ba.coarsen(refratio);  // index space is now coarsened

                // Number of guard cells to fill on coarse patch and number of components
                const amrex::IntVect ngrow = (rho[lev].nGrowVect() + refratio - 1) / refratio;  // round up int division
                const int ncomp = 1;  // rho is a scalar
                amrex::MultiFab rho_cp(ba, rho[lev].DistributionMap(), ncomp, ngrow);

                // coarsen the data
                ablastr::coarsen::average::Coarsen(
                    rho_cp,
                    rho[lev],
                    refratio
                );

                // add to the lower level
                const bool do_single_precision_comms = false;
                const amrex::Periodicity& crse_period = this->Geom(lev - 1).periodicity();
                ablastr::utils::communication::ParallelAdd(
                    rho[lev-1],
                    rho_cp,
                    0,
                    0,
                    1,
                    rho_cp.nGrowVect(),
                    rho[lev-1].nGrowVect(),
                    do_single_precision_comms,
                    crse_period
                );
            }

            // TODO: implement charge filters here
            // note: we do this after SyncRho, because the physical (dx) size of
            //       the stencil is different for each level
            // note: we might not be able to do this after SumBoundary because we would then
            //       not have valid filtered values in the guard. We access the guard
            //       when we sum contributions from particles close to the
            //       MR border between levels.

            // start async charge communication for this level
            rho_at_level.SumBoundary_nowait();
            //int const comp = 0;
            //rho_at_level.SumBoundary_nowait(comp, comp, rho_at_level.nGrowVect());

        } // lev

        // finalize communication
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            amrex::MultiFab & rho_at_level = rho.at(lev);
            rho_at_level.SumBoundary_finish();
        }
    }
} // namespace impactx
