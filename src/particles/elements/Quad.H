/* Copyright 2021 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_QUAD_H
#define IMPACTX_QUAD_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>


namespace impactx
{
    struct Quad
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A Quadrupole magnet
         *
         * @param ds Segment length in m.
         * @param k  quadrupole strength in 1/m
         */
        Quad( amrex::ParticleReal const ds, amrex::ParticleReal const k )
        : m_ds(ds), m_k(k)
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

            // advance position and momentum (drift)
            p.pos(0) = cos(m_k*m_ds)*x + sin(m_k*m_ds)/m_k*px;
            // px = px;
            px = -m_k*sin(m_k*m_ds)*x + cos(m_k*m_ds)*px;

            p.pos(1) = cosh(m_k*m_ds)*y + sinh(m_k*m_ds)/m_k*py;
            py = m_k*sinh(m_k*m_ds)*y + cosh(m_k*m_ds)*py;

            p.pos(2) = t + (m_ds/betgam2)*pt;
            // pt = pt;

        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_k; //! quadrupole strength in 1/m
    };

} // namespace impactx

#endif // IMPACTX_QUAD_H
