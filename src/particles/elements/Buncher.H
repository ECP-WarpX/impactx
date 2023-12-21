/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_BUNCHER_H
#define IMPACTX_BUNCHER_H

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
    struct Buncher
    : public elements::BeamOptic<Buncher>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "Buncher";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A short RF cavity element at zero crossing for bunching
         *
         * @param V Normalized RF voltage drop V = Emax*L/(c*Brho)
         * @param k Wavenumber of RF in 1/m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error around the s axis (in the x-y plane) [degrees]
         */
        Buncher (
            amrex::ParticleReal V,
            amrex::ParticleReal k,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_V(V), m_k(k)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a buncher functor, so that a variable of this type can be used like a
         *  buncher function.
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
                RefPart const & refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // access reference particle values to find (beta*gamma)^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum
            p.pos(RealAoS::x) = x;
            pxout = px + m_k*m_V/(2.0_prt*betgam2)*x;

            p.pos(RealAoS::y) = y;
            pyout = py + m_k*m_V/(2.0_prt*betgam2)*y;

            p.pos(RealAoS::t) = t;
            ptout = pt - m_k*m_V*t;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

    private:
        amrex::ParticleReal m_V; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_k; //! RF wavenumber in 1/m.
    };

} // namespace impactx

#endif // IMPACTX_BUNCHER_H
