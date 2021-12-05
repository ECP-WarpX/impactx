/* Copyright 2021 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DRIFT_H
#define IMPACTX_DRIFT_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct Drift
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A drift
         *
         * @param ds Segment length in m.
         */
        Drift( amrex::ParticleReal const ds )
        : m_ds(ds)
        {
        }

        /** This is a drift functor, so that a variable of this type can be used like a drift function.
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
                amrex::ParticleReal & pt) const {

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // intermediate values (universal for all elements - needs tracking for RF elements)
            amrex::ParticleReal const betgam = 2.0; // todo: make this an input parameter (reference energy -> reference betagam)
            amrex::ParticleReal const betgam2 = pow(betgam, 2);

            // advance position and momentum (drift)
            p.pos(0) = x + m_ds * px;
            // px = px;
            p.pos(1) = y + m_ds * py;
            // py = py;
            p.pos(2) = t + (m_ds/betgam2) * pt;
            // pt = pt;
        }

    private:
        amrex::ParticleReal m_ds;
    };

} // namespace impactx

#endif // IMPACTX_DRIFT_H
