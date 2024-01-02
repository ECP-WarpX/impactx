/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIPEDGE_H
#define IMPACTX_DIPEDGE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct DipEdge
    : public elements::BeamOptic<DipEdge>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "DipEdge";
        using PType = ImpactXParticleContainer::ParticleType;

        /** Edge focusing associated with bend entry or exit
         *
         * This model assumes a first-order effect of nonzero gap.
         * Here we use the linear fringe field map, given to first order in g/rc (gap / radius of curvature).
         *
         * References:
         *   K. L. Brown, SLAC Report No. 75 (1982).
         *   K. Hwang and S. Y. Lee, PRAB 18, 122401 (2015).
         *
         * @param psi Pole face angle in rad.
         * @param rc Radius of curvature in m.
         * @param g Gap parameter in m.
         * @param K2 Fringe field integral (unitless).
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        DipEdge (
            amrex::ParticleReal psi,
            amrex::ParticleReal rc,
            amrex::ParticleReal g,
            amrex::ParticleReal K2,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_psi(psi), m_rc(rc), m_g(g), m_K2(K2)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

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
                [[maybe_unused]] RefPart const & refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);

            // access reference particle values if needed

            // edge focusing matrix elements (zero gap)
            amrex::ParticleReal const R21 = tan(m_psi)/m_rc;
            amrex::ParticleReal R43 = -R21;
            amrex::ParticleReal vf = 0;

            // first-order effect of nonzero gap
            vf = (1.0_prt + pow(sin(m_psi),2))/(pow(cos(m_psi),3));
            vf *= m_g * m_K2/(pow(m_rc,2));
            R43 += vf;

            // apply edge focusing
            px = px + R21*x;
            py = py + R43*y;

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

    private:
        amrex::ParticleReal m_psi; //! pole face angle in rad
        amrex::ParticleReal m_rc; //! bend radius in m
        amrex::ParticleReal m_g; //! gap parameter in m
        amrex::ParticleReal m_K2; //! fringe field integral
    };

} // namespace impactx

#endif // IMPACTX_DIPEDGE_H
