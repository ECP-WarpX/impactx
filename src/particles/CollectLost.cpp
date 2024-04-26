/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "CollectLost.H"

#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleTransformation.H>
#include <AMReX_RandomEngine.H>


namespace impactx
{
    struct CopyAndMarkNegative
    {
        int s_index; //!< runtime index of runtime attribute in destination for position s where particle got lost
        amrex::ParticleReal s_lost; //!< position s in meters where particle got lost

        using SrcData = ImpactXParticleContainer::ParticleTileType::ConstParticleTileDataType;
        using DstData = ImpactXParticleContainer::ParticleTileType::ParticleTileDataType;

        AMREX_GPU_HOST_DEVICE
        void operator() (DstData const &dst, SrcData const &src, int src_ip, int dst_ip) const noexcept
        {
            dst.m_idcpu[dst_ip] = src.m_idcpu[src_ip];
            for (int j = 0; j < SrcData::NAR; ++j)
                dst.m_rdata[j][dst_ip] = src.m_rdata[j][src_ip];
            for (int j = 0; j < src.m_num_runtime_real; ++j)
                dst.m_runtime_rdata[j][dst_ip] = src.m_runtime_rdata[j][src_ip];

            // unused: integer compile-time or runtime attributes
            //for (int j = 0; j < SrcData::NAI; ++j)
            //    dst.m_idata[j][dst_ip] = src.m_idata[j][src_ip];
            //for (int j = 0; j < src.m_num_runtime_int; ++j)
            //    dst.m_runtime_idata[j][dst_ip] = src.m_runtime_idata[j][src_ip];

            // flip id to positive in destination
            amrex::ParticleIDWrapper{dst.m_idcpu[dst_ip]}.make_valid();

            // remember the current s of the ref particle when lost
            dst.m_runtime_rdata[s_index][dst_ip] = s_lost;
        }
    };

    void collect_lost_particles (ImpactXParticleContainer& source)
    {
        BL_PROFILE("impactX::collect_lost_particles");

        using SrcData = ImpactXParticleContainer::ParticleTileType::ConstParticleTileDataType;

        ImpactXParticleContainer& dest = *source.GetLostParticleContainer();
        const int s_runtime_index = dest.GetRealCompIndex("s_lost") - dest.NArrayReal;

        RefPart const ref_part = source.GetRefParticle();
        auto const s_lost = ref_part.s;

        // have to resize here, not in the constructor because grids have not
        // been built when constructor was called.
        dest.reserveData();
        dest.resizeData();

        // copy all particles marked with a negative ID from source to destination
        int const nLevel = source.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
            auto& plevel = source.GetParticles(lev);

            // TODO: is it safe to add OpenMP parallelism here?
            for (ParIt pti(source, lev); pti.isValid(); ++pti) {
                auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                if (plevel.find(index) == plevel.end()) continue;

                auto& ptile_source = plevel.at(index);
                auto const np = ptile_source.numParticles();
                if (np == 0) continue;  // no particles in source tile

                // we will copy particles that were marked as lost, with a negative id
                auto const predicate = [] AMREX_GPU_HOST_DEVICE (const SrcData& src, int ip)
                /* NVCC 11.3.109 chokes in C++17 on this: noexcept */
                {
                    return !amrex::ConstParticleIDWrapper{src.m_idcpu[ip]}.is_valid();
                };

                auto& ptile_dest = dest.DefineAndReturnParticleTile(
                        lev, pti.index(), pti.LocalTileIndex());

                // count how many particles we will copy
                amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
                amrex::ReduceData<int> reduce_data(reduce_op);
                {
                    auto const src_data = ptile_source.getConstParticleTileData();

                    reduce_op.eval(np, reduce_data, [=] AMREX_GPU_HOST_DEVICE (int ip)
                    {
                        return predicate(src_data, ip);
                    });
                }
                int const np_to_move = amrex::get<0>(reduce_data.value());
                if (np_to_move == 0) continue;  // no particles to move from source tile

                // allocate memory in destination
                int const dst_index = ptile_dest.numParticles();
                ptile_dest.resize(dst_index + np_to_move);

                // copy particles
                //   skipped in loop below: integer compile-time or runtime attributes
                AMREX_ALWAYS_ASSERT(SrcData::NAI == 0);
                AMREX_ALWAYS_ASSERT(ptile_source.NumRuntimeIntComps() == 0);

                //   first runtime attribute in destination is s position where particle got lost
                AMREX_ALWAYS_ASSERT(dest.NumRuntimeRealComps() > 0);

                amrex::filterAndTransformParticles(
                    ptile_dest,
                    ptile_source,
                    predicate,
                    CopyAndMarkNegative{s_runtime_index, s_lost},
                    0,
                    dst_index
                );

                // remove particles with negative ids in source
                amrex::removeInvalidParticles(ptile_source);
            } // particle tile loop
        } // lev
    }
} // namespace impactx
