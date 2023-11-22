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
#include <AMReX_ParticleTransformation.H>
#include <AMReX_RandomEngine.H>


namespace impactx
{
    void collect_lost_particles (ImpactXParticleContainer& source)
    {
        BL_PROFILE("impactX::collect_lost_particles");

        using SrcData = ImpactXParticleContainer::ParticleTileType::ConstParticleTileDataType;
        using DstData = ImpactXParticleContainer::ParticleTileType::ParticleTileDataType;

        ImpactXParticleContainer& dest = *source.GetLostParticleContainer();

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
                auto predicate = [] AMREX_GPU_HOST_DEVICE (const SrcData& src, int ip, const amrex::RandomEngine& /*engine*/) noexcept
                {
                    return src.id(ip) < 0;
                };

                auto& ptile_dest = dest.DefineAndReturnParticleTile(
                        lev, pti.index(), pti.LocalTileIndex());

                // count how many particles we will copy
                amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
                amrex::ReduceData<int> reduce_data(reduce_op);
                {
                    const auto src_data = ptile_source.getConstParticleTileData();

                    const amrex::RandomEngine rng{};  // unused
                    reduce_op.eval(np, reduce_data, [=] AMREX_GPU_HOST_DEVICE (int ip) {
                        return predicate(src_data, ip, rng) ? 1 : 0;
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

                //   first runtime attribute in destination is s position when particle got lost
                int const s_index = dest.NumRuntimeRealComps() - 1;
                auto copy_and_mark_negative = [&s_index, &s_lost](DstData& dst, const SrcData& src, int src_ip, int dst_ip) noexcept
                {
                    dst.m_aos[dst_ip] = src.m_aos[src_ip];

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
                    dst.id(dst_ip) = amrex::Math::abs(dst.id(dst_ip));

                    // remember the current s of the ref particle when lost
                    dst.m_runtime_rdata[s_index][dst_ip] = s_lost;
                };

                amrex::filterAndTransformParticles(
                    ptile_dest,
                    ptile_source,
                    predicate,
                    copy_and_mark_negative,
                    0,
                    dst_index
                );

                // remove particles with negative ids in source
                {
                    int n_removed = 0;
                    auto ptile_src_data = ptile_source.getParticleTileData();
                    for (int ip = 0; ip < np; ++ip)
                    {
                        if (ptile_source.id(ip) < 0)
                            n_removed++;
                        else
                        {
                            if (n_removed > 0)
                            {
                                // move down
                                int const new_index = ip - n_removed;

                                ptile_src_data.m_aos[new_index] = ptile_src_data.m_aos[ip];

                                for (int j = 0; j < SrcData::NAR; ++j)
                                    ptile_src_data.m_rdata[j][new_index] = ptile_src_data.m_rdata[j][ip];
                                for (int j = 0; j < ptile_src_data.m_num_runtime_real; ++j)
                                    ptile_src_data.m_runtime_rdata[j][new_index] = ptile_src_data.m_runtime_rdata[j][ip];

                                // unused: integer compile-time or runtime attributes
                                //for (int j = 0; j < SrcData::NAI; ++j)
                                //    dst.m_idata[j][new_index] = src.m_idata[j][ip];
                                //for (int j = 0; j < ptile_src_data.m_num_runtime_int; ++j)
                                //    dst.m_runtime_idata[j][new_index] = src.m_runtime_idata[j][ip];
                            }
                        }
                    }
                    AMREX_ALWAYS_ASSERT(np_to_move == n_removed);
                    ptile_source.resize(np - n_removed);
                }

            } // particle tile loop
        } // lev
    }
} // namespace impactx
