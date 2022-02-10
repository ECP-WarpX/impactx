/* Copyright 2022 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"

#include <ablastr/particles/DepositCharge.H>

#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParticleTile.H>

#include <array>


namespace impactx
{
    void
    ImpactXParticleContainer::DepositCharge (amrex::MultiFab & rho, int icomp)
    {
        // loop over refinement levels
        int const nLevel = this->maxLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            // reset the values in rho to zero
            rho.setVal(0., icomp, 1, rho.nGrowVect());

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
                    const int np = pti.numParticles();

                    // preparing access to particle data: SoA of Reals
                    auto &soa_real = pti.GetStructOfArrays().GetRealData();
                    // after https://github.com/ECP-WarpX/WarpX/pull/2838 add const:
                    auto const wp = soa_real[RealSoA::w];
                    int const * const AMREX_RESTRICT ion_lev = nullptr;

                    // mesh refinement ratio between lev and depos_lev (currently none)
                    auto const ref_ratio = amrex::IntVect(AMREX_D_DECL(1, 1, 1));
                    amrex::Box const tilebox = pti.tilebox();

                    // physical lower corner of the current box
                    amrex::RealBox const grid_box{tilebox, gm.CellSize(), gm.ProbLo()};
                    amrex::Real const * const xyzmin_ptr = grid_box.lo();
                    std::array<amrex::Real, 3> const xyzmin = {xyzmin_ptr[0], xyzmin_ptr[1], xyzmin_ptr[2]};

                    // pointer to costs data (currently unused)
                    //amrex::LayoutData<amrex::Real>* costs = WarpX::getCosts(lev);
                    //amrex::Real* cost = costs ? &((*costs)[pti.index()]) : nullptr;
                    amrex::Real * const cost = nullptr;
                    long const load_balance_costs_update_algo = 0;

                    // used for MR when we want to deposit for a subset of the particles on the level in the
                    // current box; with offset, we start at a later particle index
                    long const offset = 0;
                    long const np_to_depose = np;
                    int const depos_lev = lev;

                    amrex::ParticleReal const charge = 1.0; // TODO once we implement charge

                    // particle shapes and needed deposition guards
                    auto const[nox, noy, noz] = std::array<int, 3>{m_particle_shape, m_particle_shape, m_particle_shape};
                    auto const ng_rho = amrex::IntVect{nox + 1, noy + 1, noz + 1};

                    // cell size of the mesh to deposit to
                    std::array<amrex::Real, 3> const &dx = {gm.CellSize(0), gm.CellSize(1), gm.CellSize(2)};

                    // number of components to deposit
                    int const nc = 1;

                    // RZ modes (unused)
                    int const n_rz_azimuthal_modes = 0;

                    // call amrex::Gpu::synchronize() for tiny profiler regions
                    bool const do_device_synchronize = true;

                    ablastr::particles::deposit_charge<ImpactXParticleContainer>
                            (pti, wp, ion_lev, &rho, icomp, nc, offset, np_to_depose,
                             local_rho_fab, lev, depos_lev, charge,
                             nox, noy, noz, ng_rho, dx, xyzmin, ref_ratio,
                             cost, n_rz_azimuthal_modes,
                             load_balance_costs_update_algo,
                             do_device_synchronize);
                }
            }
        }
    }
} // namespace impactx
