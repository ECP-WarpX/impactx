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


namespace impactx
{
    void collect_lost_particles (ImpactXParticleContainer& source, ImpactXParticleContainer& dest)
    {
        using SrcData = ImpactXParticleContainer::ParticleTileType::ConstParticleTileDataType;

        // copy all particles marked with a negative ID from source to destination
        bool const local = true;
        dest.copyParticles(
                source,
                [=] AMREX_GPU_HOST_DEVICE(const SrcData &src, int ip) {
                    return src.id(ip) < 0;
                },
                local
        );

        // flip IDs back to positive in destination
        int const nLevel = dest.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(dest, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto &aos = pti.GetArrayOfStructs();
                PType *AMREX_RESTRICT aos_ptr = aos().dataPtr();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {
                    PType &p = aos_ptr[i];
                    p.id() = -p.id();
                });
            }
        }

        // TODO remove particles with negative ids in source
    }
} // namespace impactx
