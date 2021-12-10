/* Copyright 2021 Axel Huebl, Remi Lehe
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"

#include <ablastr/particles/DepositCharge.H>
#include <ablastr/particles/ParticleMoments.H>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleTile.H>

#include <array>


namespace impactx
{
    ImpactXParticleContainer::ImpactXParticleContainer (amrex::AmrCore* amr_core)
        : amrex::ParticleContainer<0, 0, RealSoA::nattribs, IntSoA::nattribs>(amr_core->GetParGDB())
    {
        SetParticleSize();

        // particle shapes
        amrex::ParmParse pp_algo("algo");
        pp_algo.get("particle_shape", m_particle_shape);
        if (m_particle_shape < 1 || m_particle_shape > 3)
        {
            amrex::Abort("algo.particle_shape can be only 1, 2, or 3");
        }
    }

    void
    ImpactXParticleContainer::AddNParticles (int lev,
                                             amrex::Vector<amrex::ParticleReal> const & x,
                                             amrex::Vector<amrex::ParticleReal> const & y,
                                             amrex::Vector<amrex::ParticleReal> const & z,
                                             amrex::Vector<amrex::ParticleReal> const & px,
                                             amrex::Vector<amrex::ParticleReal> const & py,
                                             amrex::Vector<amrex::ParticleReal> const & pz)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0, "AddNParticles: only lev=0 is supported yet.");
        AMREX_ALWAYS_ASSERT(x.size() == y.size());
        AMREX_ALWAYS_ASSERT(x.size() == z.size());

        // number of particles to add
        int const np = x.size();

        // have to resize here, not in the constructor because grids have not
        // been built when constructor was called.
        reserveData();
        resizeData();

        auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);

        /* Create a temporary tile to obtain data from simulation. This data
         * is then copied to the permanent tile which is stored on the particle
         * (particle_tile).
         */
        using PinnedTile = amrex::ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt,
                amrex::PinnedArenaAllocator>;
        PinnedTile pinned_tile;
        pinned_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());

        for (int i = 0; i < np; i++)
        {
            ParticleType p;
            p.id() = ParticleType::NextID();
            p.cpu() = amrex::ParallelDescriptor::MyProc();
            p.pos(0) = x[i];
            p.pos(1) = y[i];
            p.pos(2) = z[i];
            // write position, creating cpu id, and particle id
            pinned_tile.push_back(p);
        }

        // write Real attributes (SoA) to particle initialized zero
        DefineAndReturnParticleTile(0, 0, 0);

        pinned_tile.push_back_real(RealSoA::ux, np, 0.0);
        pinned_tile.push_back_real(RealSoA::uy, np, 0.0);
        pinned_tile.push_back_real(RealSoA::pt, np, 0.0);
        pinned_tile.push_back_real(RealSoA::ux, *px.cbegin(), *px.cend());
        pinned_tile.push_back_real(RealSoA::uy, *py.cbegin(), *py.cend());
        pinned_tile.push_back_real(RealSoA::pt, *pz.cbegin(), *pz.cend());

        //the following should be updated
        pinned_tile.push_back_real(RealSoA::t, np, 0.0);
        pinned_tile.push_back_real(RealSoA::q_m, np, 0.0);
        pinned_tile.push_back_real(RealSoA::w, np, 1.0/np);

        /* Redistributes particles to their respective tiles (spatial bucket
         * sort per box over MPI ranks)
         */
        auto old_np = particle_tile.numParticles();
        auto new_np = old_np + pinned_tile.numParticles();
        particle_tile.resize(new_np);
        amrex::copyParticles(
                particle_tile, pinned_tile, 0, old_np, pinned_tile.numParticles());
//        Redistribute(); // TODO

    }

    void
    ImpactXParticleContainer::SetRefParticle (RefPart const refpart)
    {
        m_refpart = refpart;
    }

    RefPart
    ImpactXParticleContainer::GetRefParticle () const
    {
        return m_refpart;
    }

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MinAndMaxPositions ()
    {
        return ablastr::particles::MinAndMaxPositions(*this);
    }

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MeanAndStdPositions ()
    {
        return ablastr::particles::MeanAndStdPositions<
            ImpactXParticleContainer, RealSoA::w
        >(*this);
    }

    void
    DepositCharge (ImpactXParticleContainer & pc, amrex::MultiFab* rho, int icomp) // ImpactXParticleContainer & pc)
    {
        // loop over refinement levels
        int const nLevel = pc.maxLevel();
        for (int lev = 0; lev < nLevel; ++lev) {
            // get simulation geometry information
            //const amrex::Geometry& gm = this->Geom(lev);
            //const auto prob_lo = gm.ProbLo();

            amrex::FArrayBox local_rho_fab;

            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                //const auto t_lev = pti.GetLevel();
                //const auto index = pti.GetPairIndex();
                // ...

                // preparing access to particle data: SoA of Reals
                auto &soa_real = pti.GetStructOfArrays().GetRealData();
                //amrex::ParticleReal *const AMREX_RESTRICT wp = soa_real[RealSoA::w];
                auto wp = soa_real[RealSoA::w];
                int const * const AMREX_RESTRICT ion_lev = nullptr;

                // DepositCharge(pti, wp, ion_lev, rho[lev].get(), icomp, 0, np, thread_num, lev, lev);

                amrex::IntVect ref_ratio = amrex::IntVect(AMREX_D_DECL(1, 1, 1));
                amrex::Box tilebox = pti.tilebox();

                // TODO: physical lower corner of the current box
                amrex::Geometry const &gm = pc.Geom(lev);
                amrex::RealBox const grid_box{tilebox, gm.CellSize(), gm.ProbLo()};
                amrex::Real const * const xyzmin_ptr = grid_box.lo();
                std::array<amrex::Real, 3> const xyzmin = {xyzmin_ptr[0], xyzmin_ptr[1], xyzmin_ptr[2]};

                // pointer to costs data
                //amrex::LayoutData<amrex::Real>* costs = WarpX::getCosts(lev);
                //amrex::Real* cost = costs ? &((*costs)[pti.index()]) : nullptr;
                amrex::Real *cost = nullptr;

                // used for MR when we want to deposit for a subset of the particles on the level in the
                // current box; with offset, we start at a later particle index
                const long offset = 0;
                const long np_to_depose = np;
                int depos_lev = lev;

                amrex::ParticleReal const charge = 1.0; // fixme

                auto const[nox, noy, noz] = std::array<int, 3>{2, 2, 2};
                auto const ng_rho = amrex::IntVect{nox + 1, noy + 1, noz + 1};

                const std::array<amrex::Real, 3> &dx = {1., 1., 1.};

                // number of components to deposit
                const int nc = 1;

                const int n_rz_azimuthal_modes = 0;
                const long load_balance_costs_update_algo = 0;

                // call amrex::Gpu::synchronize() for tiny profiler regions
                const bool do_device_synchronize = true;

                ablastr::particles::deposit_charge<ImpactXParticleContainer>
                        (pti, wp, ion_lev, rho, icomp, nc, offset, np_to_depose,
                         local_rho_fab, lev, depos_lev, charge,
                         nox, noy, noz, ng_rho, dx, xyzmin, ref_ratio,
                         cost, n_rz_azimuthal_modes,
                         load_balance_costs_update_algo,
                         do_device_synchronize);
            }
        }
    }
} // namespace impactx
