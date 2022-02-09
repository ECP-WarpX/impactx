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

#include <cmath>


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

        /** This is a quad functor, so that a variable of this type can be used like a quad function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                PType& p,
                amrex::ParticleReal & px,
                amrex::ParticleReal & py,
                amrex::ParticleReal & pt,
                RefPart const refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

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
