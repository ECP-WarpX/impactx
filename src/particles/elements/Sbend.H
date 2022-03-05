/* Copyright 2021 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SBEND_H
#define IMPACTX_SBEND_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


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
        Sbend( amrex::ParticleReal const ds, amrex::ParticleReal const rc)
        : m_ds(ds), m_rc(rc)
        {
        }

        /** This is a sbend functor, so that a variable of this type can be used like a sbend function.
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

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;
            amrex::ParticleReal const bet = sqrt(betgam2/(1.0_prt + betgam2));
            amrex::ParticleReal const theta = m_ds/m_rc;

            // advance position and momentum (sector bend)

            p.pos(0) = cos(theta)*x + m_rc*sin(theta)*px
                     - (m_rc/bet)*(1.0_prt - cos(theta))*pt;

            pxout = -sin(theta)/m_rc*x + cos(theta)*px - sin(theta)/bet*pt;

            p.pos(1) = y + m_rc*theta*py;

            // pyout = py;

            p.pos(2) = sin(theta)/bet*x + m_rc/bet*(1.0_prt - cos(theta))*px + t
                     + m_rc*(-theta+sin(theta)/(bet*bet))*pt;

            // ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_rc; //! bend radius in m
    };

} // namespace impactx

#endif // IMPACTX_SBEND_H
