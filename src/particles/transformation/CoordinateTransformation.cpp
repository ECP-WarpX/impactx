/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "CoordinateTransformation.H"

#include "T2Z.H"
#include "Z2T.H"

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal

#include <cmath>


namespace impactx
{
namespace transformation {
    void CoordinateTransformation (ImpactXParticleContainer &pc,
                                   Direction const &direction) {
        using namespace amrex::literals; // for _rt and _prt

        // preparing to access reference particle data: RefPart
        RefPart ref_part = pc.GetRefParticle();
        amrex::ParticleReal const pd = ref_part.pt;  // Design value of pt/mc2 = -gamma

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto &aos = pti.GetArrayOfStructs();
                PType *AMREX_RESTRICT aos_ptr = aos().dataPtr();

                // preparing access to particle data: SoA of Reals
                auto &soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal *const AMREX_RESTRICT part_px = soa_real[RealSoA::ux].dataPtr();
                amrex::ParticleReal *const AMREX_RESTRICT part_py = soa_real[RealSoA::uy].dataPtr();
                amrex::ParticleReal *const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

                if( direction == Direction::T2Z) {
                    // Design value of pz/mc = beta*gamma
                    amrex::ParticleReal const pzd = sqrt(pow(pd, 2) - 1.0);

                    T2Z t2z(pzd);
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {
                        // access AoS data such as positions and cpu/id
                        PType &p = aos_ptr[i];

                        // access SoA Real data
                        amrex::ParticleReal &px = part_px[i];
                        amrex::ParticleReal &py = part_py[i];
                        amrex::ParticleReal &pt = part_pt[i];

                        t2z(p, px, py, pt);
                    });
                } else {
                    amrex::ParticleReal const ptd = pd;  // Design value of pt/mc2 = -gamma.
                    Z2T z2t(ptd);
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {
                        // access AoS data such as positions and cpu/id
                        PType &p = aos_ptr[i];

                        // access SoA Real data
                        amrex::ParticleReal &px = part_px[i];
                        amrex::ParticleReal &py = part_py[i];
                        amrex::ParticleReal &pt = part_pt[i];

                        z2t(p, px, py, pt);
                    });
                }
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }
} // namespace transformation
} // namespace impactx
