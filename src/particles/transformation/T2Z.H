/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_T2Z_H
#define IMPACTX_T2Z_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
namespace transformation
{
    struct T2Z
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A T-to-Z transformation
         *
         * @param pzd Design value of pz/mc = beta*gamma.
         */
        T2Z( amrex::ParticleReal const pzd )
        : m_pzd(pzd)
        {
        }

        /** This is a t-to-z map, so that a variable of this type can be used like a t-to-z function.
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

            // arguments for ptd and pt
            amrex::ParticleReal const argd = 1.0_prt + pow(m_pzd, 2);
            amrex::ParticleReal const arg = 1.0_prt + pow(m_pzd+pt, 2) + pow(px, 2) + pow(py, 2);

            // check for validity of arguments
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(argd > 0.0_prt, "invalid ptd arg (<=0)");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(arg > 0.0_prt, "invalid pt arg (<=0)");

            // transform position and momentum (t2z)
            //   (universal for all elements - needs tracking for RF elements)
            amrex::ParticleReal const ptdf = argd > 0.0_prt ? -sqrt(argd) : -1.0_prt;
            amrex::ParticleReal const ptf = arg > 0.0_prt ? -sqrt(arg) : -1.0_prt;

            p.pos(0) = x - px*t/(m_pzd+pt);
            // px = px;
            p.pos(1) = y - py*t/(m_pzd+pt);
            // py = py;
            p.pos(2) = ptf*t/(m_pzd+pt);  // This now represents t.
            pt = ptf - ptdf;                    // This now represents pt.
        }

    private:
        amrex::ParticleReal m_pzd;
    };

} // namespace transformation
} // namespace impactx

#endif // IMPACTX_T2Z_H
