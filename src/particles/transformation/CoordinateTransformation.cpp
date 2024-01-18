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


namespace impactx::transformation
{
    void CoordinateTransformation (ImpactXParticleContainer & pc,
                                   CoordSystem direction)
    {
        BL_PROFILE("impactx::transformation::CoordinateTransformation");
        using namespace amrex::literals; // for _rt and _prt

        if (direction == CoordSystem::s) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pc.GetCoordSystem() == CoordSystem::t, "Already in fixed s coordinates!");
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pc.GetCoordSystem() == CoordSystem::s, "Already in fixed t coordinates!");
        }

        // preparing to access reference particle data: RefPart
        RefPart const ref_part = pc.GetRefParticle();
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

                // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto &aos = pti.GetArrayOfStructs();
                PType *AMREX_RESTRICT aos_ptr = aos().dataPtr();

                // preparing access to particle data: SoA of Reals
                auto &soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal *const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal *const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();

                if (direction == CoordSystem::s) {
                    BL_PROFILE("impactx::transformation::CoordinateTransformation::to_fixed_s");

                    amrex::ParticleReal *const AMREX_RESTRICT part_pz = soa_real[RealSoA::pz].dataPtr();

                    // Design value of pz/mc = beta*gamma
                    amrex::ParticleReal const pzd = sqrt(pow(pd, 2) - 1.0);

                    ToFixedS const to_s(pzd);
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {
                        // access AoS data such as positions and cpu/id
                        PType &p = aos_ptr[i];

                        // access SoA Real data
                        amrex::ParticleReal &px = part_px[i];
                        amrex::ParticleReal &py = part_py[i];
                        amrex::ParticleReal &pz = part_pz[i];

                        to_s(p, px, py, pz);
                    });
                } else {
                    BL_PROFILE("impactx::transformation::CoordinateTransformation::to_fixed_t");

                    amrex::ParticleReal *const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

                    amrex::ParticleReal const ptd = pd;  // Design value of pt/mc2 = -gamma.
                    ToFixedT const to_t(ptd);
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {
                        // access AoS data such as positions and cpu/id
                        PType &p = aos_ptr[i];

                        // access SoA Real data
                        amrex::ParticleReal &px = part_px[i];
                        amrex::ParticleReal &py = part_py[i];
                        amrex::ParticleReal &pt = part_pt[i];

                        to_t(p, px, py, pt);
                    });
                }
            } // end loop over all particle boxes
        } // env mesh-refinement level loop

        // update coordinate system meta data
        pc.SetCoordSystem(direction);
    }
} // namespace impactx::transformation
