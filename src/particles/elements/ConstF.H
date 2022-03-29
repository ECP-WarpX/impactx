/* Copyright 2022 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CONSTF_H
#define IMPACTX_CONSTF_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct ConstF
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A linear Constant Focusing element
         *
         * @param ds Segment length in m.
         * @param kx Focusing strength for x in 1/m.
         * @param ky Focusing strength for y in 1/m.
         * @param kt Focusing strength for t in 1/m.
         *
         */
        ConstF( amrex::ParticleReal const ds, amrex::ParticleReal const kx,
                amrex::ParticleReal const ky, amrex::ParticleReal const kt )
        : m_ds(ds), m_kx(kx), m_ky(ky), m_kt(kt)
        {
        }

        /** This is a constf functor, so that a variable of this type can be used like a
         *  constf function.
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
                [[maybe_unused]] RefPart const refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum
            p.pos(0) = cos(m_kx*m_ds)*x + sin(m_kx*m_ds)/m_kx*px;
            pxout = -m_kx*sin(m_kx*m_ds)*x + cos(m_kx*m_ds)*px;

            p.pos(1) = cos(m_ky*m_ds)*y + sin(m_ky*m_ds)/m_ky*py;
            pyout = -m_ky*sin(m_ky*m_ds)*y + cos(m_ky*m_ds)*py;

            p.pos(2) = cos(m_kt*m_ds)*t + sin(m_kt*m_ds)/(betgam2*m_kt)*pt;
            ptout = -(m_kt*betgam2)*sin(m_kt*m_ds)*t + cos(m_kt*m_ds)*pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_kx; //! focusing x strength in 1/m
        amrex::ParticleReal m_ky; //! focusing y strength in 1/m
        amrex::ParticleReal m_kt; //! focusing t strength in 1/m

    };

} // namespace impactx

#endif // IMPACTX_CONSTF_H
