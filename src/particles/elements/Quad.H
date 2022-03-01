/* Copyright 2021 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_QUAD_H
#define IMPACTX_QUAD_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
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
                PType& AMREX_RESTRICT p,
                amrex::ParticleReal & AMREX_RESTRICT px,
                amrex::ParticleReal & AMREX_RESTRICT py,
                amrex::ParticleReal & AMREX_RESTRICT pt,
                RefPart const refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            amrex::ParticleReal const omega = sqrt(abs(m_k));

            if(m_k > 0.0) {
               // advance position and momentum (focusing quad)
               p.pos(0) = cos(omega*m_ds)*x + sin(omega*m_ds)/omega*px;
               px = -omega*sin(omega*m_ds)*x + cos(omega*m_ds)*px;

               p.pos(1) = cosh(omega*m_ds)*y + sinh(omega*m_ds)/omega*py;
               py = omega*sinh(omega*m_ds)*y + cosh(omega*m_ds)*py;

               p.pos(2) = t + (m_ds/betgam2)*pt;
               // pt = pt;
            } else {
               // advance position and momentum (defocusing quad)
               p.pos(0) = cosh(omega*m_ds)*x + sinh(omega*m_ds)/omega*px;
               px = omega*sinh(omega*m_ds)*x + cosh(omega*m_ds)*px;

               p.pos(1) = cos(omega*m_ds)*y + sin(omega*m_ds)/omega*py;
               py = -omega*sin(omega*m_ds)*y + cos(omega*m_ds)*py;

               p.pos(2) = t + (m_ds/betgam2)*pt;
               // pt = pt;
            }
        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_k; //! quadrupole strength in 1/m
    };

} // namespace impactx

#endif // IMPACTX_QUAD_H
