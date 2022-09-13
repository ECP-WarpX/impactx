/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "CoordinateTransformation.H"

#include "ToFixedS.H"
#include "ToFixedT.H"

#include <AMReX_BLProfiler.H> // for BL_PROFILE
#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal

#include <cmath>


namespace impactx
{
namespace transformation {
    void CoordinateTransformation (ImpactXParticleContainer &pc,
                                   Direction const &direction)
   {
        BL_PROFILE("impactx::transformation::CoordinateTransformation");
        using namespace amrex::literals; // for _rt and _prt

        // preparing to access reference particle data: RefPart
        RefPart ref_part = pc.GetRefParticle();
        amrex::ParticleReal const pd = ref_part.pt;  // Design value of pt/mc2 = -gamma

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                using PType = ImpactXParticleContainer::ParticleType;

                // preparing access to particle data: SoA of Reals
                auto &soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal *const AMREX_RESTRICT part_px = soa_real[RealSoA::ux].dataPtr();
                amrex::ParticleReal *const AMREX_RESTRICT part_py = soa_real[RealSoA::uy].dataPtr();
                amrex::ParticleReal *const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

                if( direction == Direction::to_fixed_s) {
                    BL_PROFILE("impactx::transformation::CoordinateTransformation::to_fixed_s");
                    // Design value of pz/mc = beta*gamma
                    amrex::ParticleReal const pzd = sqrt(pow(pd, 2) - 1.0);

                    ToFixedS const to_s(pzd);
                    auto tile_data = pti.GetParticleTile().getParticleTileData();

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {

                        PType p(tile_data,i);

                        // access SoA Real data
                        amrex::ParticleReal &px = part_px[i];
                        amrex::ParticleReal &py = part_py[i];
                        amrex::ParticleReal &pz = part_pt[i];

                        to_s(p, px, py, pz);
                    });
                } else {
                    BL_PROFILE("impactx::transformation::CoordinateTransformation::to_fixed_t");
                    amrex::ParticleReal const ptd = pd;  // Design value of pt/mc2 = -gamma.
                    ToFixedT const to_t(ptd);
                    auto tile_data = pti.GetParticleTile().getParticleTileData();
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {
                        // access AoS data such as positions and cpu/id
                        PType p(tile_data,i);

                        // access SoA Real data
                        amrex::ParticleReal &px = part_px[i];
                        amrex::ParticleReal &py = part_py[i];
                        amrex::ParticleReal &pt = part_pt[i];

                        to_t(p, px, py, pt);
                    });
                }
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }
} // namespace transformation
} // namespace impactx
