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

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MinAndMaxPositions ()
    {
        using PType = ImpactXParticleContainer::SuperParticleType;

        // Get min and max for the local rank
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

        // Get min and max across all ranks
        std::vector<amrex::ParticleReal> xyz_min = {
            amrex::get<0>(r),
            amrex::get<1>(r),
            amrex::get<2>(r)
        };
        amrex::ParallelDescriptor::ReduceRealMin(xyz_min.data(), xyz_min.size());
        std::vector<amrex::ParticleReal> xyz_max = {
            amrex::get<3>(r),
            amrex::get<4>(r),
            amrex::get<5>(r)
        };
        amrex::ParallelDescriptor::ReduceRealMax(xyz_max.data(), xyz_max.size());

        return {xyz_min[0], xyz_min[1], xyz_min[2], xyz_max[0], xyz_max[1], xyz_max[2]};
    }

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MeanAndStdPositions ()
    {
        using PType = ImpactXParticleContainer::SuperParticleType;

        amrex::ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_ops;
        auto r = amrex::ParticleReduce<amrex::ReduceData<ParticleReal, ParticleReal, ParticleReal, ParticleReal, ParticleReal, ParticleReal, ParticleReal>>(
            *this,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<ParticleReal, ParticleReal, ParticleReal, ParticleReal, ParticleReal, ParticleReal, ParticleReal>
            {
                amrex::ParticleReal x = p.pos(0);
                amrex::ParticleReal y = p.pos(1);
                amrex::ParticleReal z = p.pos(2);
                amrex::ParticleReal w = p.rdata(RealSoA::w);

                return {x, x*x, y, y*y, z, z*z, w};
            },
            reduce_ops);

        // Reduce across MPI ranks
        std::vector<amrex::ParticleReal> data_vector = {
            amrex::get<0>(r),
            amrex::get<1>(r),
            amrex::get<2>(r),
            amrex::get<3>(r);
            amrex::get<4>(r);
            amrex::get<5>(r);
            amrex::get<6>(r);
        };
        ParallelDescriptor::ReduceRealSum(data_vector.data(), data_vector.size());

        amrex::ParticleReal w_sum = data_vector[6];
        amrex::ParticleReal x_mean = data_vector[0]/w_sum;
        amrex::ParticleReal x_std = data_vector[1]/w_sum- x_mean*x_mean;
        amrex::ParticleReal y_mean = data_vector[2]/w_sum;
        amrex::ParticleReal y_std = data_vector[3]/w_sum- x_mean*x_mean;
        amrex::ParticleReal z_mean = data_vector[4]/w_sum;
        amrex::ParticleReal z_std = data_vector[5]/w_sum- x_mean*x_mean;

        return {x_mean, x_std, y_mean, y_std, z_mean, z_std};
    }

} // namespace impactx
