/* Copyright 2021 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"

#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParticleTile.H>


namespace impactx
{
    ImpactXParticleContainer::ImpactXParticleContainer (amrex::AmrCore* amr_core)
        : amrex::ParticleContainer<0, 0, RealSoA::nattribs, IntSoA::nattribs>(amr_core->GetParGDB())
    {
       SetParticleSize();
    }

    void
    ImpactXParticleContainer::AddNParticles (int lev,
                                             amrex::Vector<amrex::ParticleReal> const & x,
                                             amrex::Vector<amrex::ParticleReal> const & y,
                                             amrex::Vector<amrex::ParticleReal> const & z)
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
        pinned_tile.push_back_real(RealSoA::t, np, 0.0);
        pinned_tile.push_back_real(RealSoA::pt, np, 0.0);
        pinned_tile.push_back_real(RealSoA::q_m, np, 0.0);
        pinned_tile.push_back_real(RealSoA::w, np, 0.0);

        /* Redistributes particles to their respective tiles (spatial bucket
         * sort per box over MPI ranks)
         */
        auto old_np = particle_tile.numParticles();
        auto new_np = old_np + pinned_tile.numParticles();
        particle_tile.resize(new_np);
        amrex::copyParticles(
                particle_tile, pinned_tile, 0, old_np, pinned_tile.numParticles());
//        Redistribute();

    }

    void
    ImpactXParticleContainer::MinAndMaxPositions (
        amrex::ParticleReal& x_min, amrex::ParticleReal& x_max,
        amrex::ParticleReal& y_min, amrex::ParticleReal& y_max,
        amrex::ParticleReal& z_min, amrex::ParticleReal& z_max )
    {
        using PType = ImpactXParticleContainer::SuperParticleType;

        amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_ops;
        auto r = amrex::ParticleReduce<amrex::ReduceData<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal>>(
            *this,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal>
            {
                amrex::ParticleReal x = p.pos(0);
                amrex::ParticleReal y = p.pos(1);
                amrex::ParticleReal z = p.pos(2);
                return {x, y, z, x, y, z};
            },
            reduce_ops);

        x_min = amrex::get<0>(r);
        y_min = amrex::get<1>(r);
        z_min = amrex::get<2>(r);
        x_max = amrex::get<3>(r);
        y_max = amrex::get<4>(r);
        z_max = amrex::get<5>(r);

        amrex::ParallelDescriptor::ReduceRealMin(x_min);
        amrex::ParallelDescriptor::ReduceRealMin(y_min);
        amrex::ParallelDescriptor::ReduceRealMin(z_min);
        amrex::ParallelDescriptor::ReduceRealMax(x_max);
        amrex::ParallelDescriptor::ReduceRealMax(y_max);
        amrex::ParallelDescriptor::ReduceRealMax(z_max);
    }

} // namespace impactx
