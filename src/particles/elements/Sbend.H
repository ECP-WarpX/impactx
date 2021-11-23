/* Copyright 2021 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SBEND_H
#define IMPACTX_SBEND_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>


namespace impactx
{
    struct Sbend
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** An ideal sector bend
         *
         * @param ds Segment length in m.
         * @param rc Radius of curvature in m.
         */
        Sbend( amrex::ParticleReal const ds, amrex::ParticleReal const rc )
        : m_ds(ds), m_rc(rc)
        {
        }

        /** ...
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
                amrex::ParticleReal & pt) {

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // intermediate values (universal for all elements - needs tracking for RF elements)
            amrex::ParticleReal const betgam = 2.0; // todo: make this an input parameter (reference energy -> reference betagam)
            amrex::ParticleReal const betgam2 = pow(betgam, 2);
            amrex::ParticleReal const bet = sqrt(betgam2/(1.0+betgam2));
            amrex::ParticleReal const theta = m_ds/m_rc;

            // advance position and momentum (sector bend)
            p.pos(0) = cos(theta)*x + m_rc*sin(theta)*px - (m_rc/bet)*(1.0-cos(theta))*pt;
            px = -sin(theta)/m_rc*x + cos(theta)*px - sin(theta)/bet*pt;

            p.pos(1) = y + m_rc*theta*py;
            // py = py;

            p.pos(2) = sin(theta)/bet*x + m_rc/bet*(1.0-cos(theta))*px + t
                     + m_rc*(-theta+sin(theta)/(bet*bet))*pt;
            // pt = pt;
        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_rc; //! bend radius in m
    };

} // namespace impactx

#endif // IMPACTX_SBEND_H
