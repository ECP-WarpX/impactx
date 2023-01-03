/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TO_FIXED_T_H
#define IMPACTX_TO_FIXED_T_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
namespace transformation
{
    struct ToFixedT
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** Transformation of particles from fixed location s to fixed time t.
         *
         * At fixed t, each particle is represented by phase space
         *  coordinates (x,y,z,px,py,pz).
         *
         * @param ptd Design value of pt/mc2 = -gamma.
         */
        ToFixedT( amrex::ParticleReal const ptd )
        : m_ptd(ptd)
        {
        }

        /** This is a s-to-t map, so that a variable of this type can be used like a
         *  s-to-t function.
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

            // compute value of reference pzd = beta*gamma
            amrex::ParticleReal const argd = -1.0_prt + pow(m_ptd, 2);
            AMREX_ASSERT_WITH_MESSAGE(argd > 0.0_prt, "invalid pzd arg (<=0)");
            amrex::ParticleReal const pzd = argd > 0.0_prt ? sqrt(argd) : 0.0_prt;

            // transform momenta to dynamic units (eg, so that momenta are
            // normalized by mc):
            px = px*pzd;
            py = py*pzd;
            pt = pt*pzd;

            // compute value of particle pz = beta*gamma
            amrex::ParticleReal const arg = -1.0_prt + pow(m_ptd+pt, 2) - pow(px, 2) - pow(py, 2);
            AMREX_ASSERT_WITH_MESSAGE(arg > 0.0_prt, "invalid pz arg (<=0)");
            amrex::ParticleReal const pz = arg > 0.0_prt ? sqrt(arg) : 0.0_prt;

            // transform position and momentum (from fixed s to fixed t)

            p.pos(0) = x + px*t/(m_ptd+pt);
            // px = px;
            p.pos(1) = y + py*t/(m_ptd+pt);
            // py = py;
            p.pos(2) = pz*t/(m_ptd+pt);  // This now represents z.
            pt = pz - pzd;                     // This now represents pz.

            // transform momenta to static units (eg, so that momenta are
            // normalized by pzd):
            px = px/pzd;
            py = py/pzd;
            pt = pt/pzd;

        }

    private:
        amrex::ParticleReal m_ptd;
    };

} // namespace transformation
} // namespace impactx

#endif // IMPACTX_TO_FIXED_T_H
