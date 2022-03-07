/* Copyright 2021-2022 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIPEDGE_H
#define IMPACTX_DIPEDGE_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct DipEdge
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** Edge focusing associated with bend entry or exit
         *
         * @param psi Pole face angle in rad.
         * @param rc Radius of curvature in m.
         * @param g Gap parameter in m.
         * @param K2 Fringe field integral (unitless).
         */
        DipEdge( amrex::ParticleReal const psi, amrex::ParticleReal const rc, amrex::ParticleReal const
        g, amrex::ParticleReal const K2)
        : m_psi(psi), m_rc(rc), m_g(g), m_K2(K2)
        {
        }

        /** This is a dipedge functor, so that a variable of this type can be used like a
         *  dipedge function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unchanged)
         * @param refpart reference particle (unused)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                PType& AMREX_RESTRICT p,
                amrex::ParticleReal & AMREX_RESTRICT px,
                amrex::ParticleReal & AMREX_RESTRICT py,
                [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT pt,
                [[maybe_unused]] RefPart const refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);

            // access reference particle values if needed

            // edge focusing matrix elements (zero gap)
            amrex::ParticleReal R21 = tan(m_psi)/m_rc;
            amrex::ParticleReal R43 = -tan(m_psi)/m_rc;
            amrex::ParticleReal vf = 0;

            // first-order effect of nonzero gap
            vf = (1.0_prt + pow(sin(m_psi),2))/(pow(cos(m_psi),3));
            vf = vf * m_g * m_K2/(pow(m_rc,2));
            R43 = R43 + vf;

            // apply edge focusing
            px = px + R21*x;
            py = py + R43*y;
        }

    private:
        amrex::ParticleReal m_psi; //! pole face angle in rad
        amrex::ParticleReal m_rc; //! bend radius in m
        amrex::ParticleReal m_g; //! gap parameter in m
        amrex::ParticleReal m_K2; //! fringe field integral
    };

} // namespace impactx

#endif // IMPACTX_DIPEDGE_H
