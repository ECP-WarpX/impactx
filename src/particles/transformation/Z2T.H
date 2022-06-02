/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_Z2T_H
#define IMPACTX_Z2T_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
namespace transformation
{
    struct Z2T
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A Z-to-T transformation
         *
         * @param ptd Design value of pt/mc2 = -gamma.
         */
        Z2T( amrex::ParticleReal const ptd )
        : m_ptd(ptd)
        {
        }

        /** This is a z-to-t map, so that a variable of this type can be used like a z-to-t function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                PType& p,
                amrex::ParticleReal & px,
                amrex::ParticleReal & py,
                amrex::ParticleReal & pt) const
        {
            using namespace amrex::literals;

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // arguments for pzd and pz
            amrex::ParticleReal const argd = -1.0_prt + pow(m_ptd, 2);
            amrex::ParticleReal const arg = -1.0_prt + pow(m_ptd+pt, 2) - pow(px, 2) - pow(py, 2);

            // check for validity of arguments
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(argd > 0.0_prt, "invalid pzd arg (<=0)");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(arg > 0.0_prt, "invalid pz arg (<=0)");

            // transform position and momentum (z2t)
            //   (universal for all elements - needs tracking for RF elements)
            amrex::ParticleReal const pzd = argd > 0.0_prt ? sqrt(argd) : 0.0_prt;
            amrex::ParticleReal const pz = arg > 0.0_prt ? sqrt(arg) : 0.0_prt;

            p.pos(0) = x + px*t/(m_ptd+pt);
            // px = px;
            p.pos(1) = y + py*t/(m_ptd+pt);
            // py = py;
            p.pos(2) = pz*t/(m_ptd+pt);  // This now represents z.
            pt = pz - pzd;                     // This now represents pz.
        }

    private:
        amrex::ParticleReal m_ptd;
    };

} // namespace transformation
} // namespace impactx

#endif // IMPACTX_Z2T_H
