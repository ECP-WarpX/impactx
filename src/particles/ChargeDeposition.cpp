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
        // loop over refinement levels
        int const nLevel = this->finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            amrex::MultiFab & rho_at_level = rho.at(lev);
            // reset the values in rho to zero
            rho_at_level.setVal(0.);

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

            // start async charge communication for this level
            rho_at_level.SumBoundary_nowait();
            //int const comp = 0;
            //rho_at_level.SumBoundary_nowait(comp, comp, rho_at_level.nGrowVect());
        }

        // finalize communication
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            amrex::MultiFab & rho_at_level = rho.at(lev);
            rho_at_level.SumBoundary_finish();
        }
    }
} // namespace impactx
