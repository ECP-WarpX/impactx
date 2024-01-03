/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SHORTRF_H
#define IMPACTX_SHORTRF_H

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
    struct ShortRF
    : public elements::BeamOptic<ShortRF>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "ShortRF";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A short RF cavity element
         *
         * @param V Normalized RF voltage V = maximum energy gain/(m*c^2)
         * @param freq RF frequency in Hz
         * @param phase RF synchronous phase of the reference particle in deg
         *     phi = 0:  maximum energy gain (on-crest)
         *     phi = -90 deg:  zero-crossing for bunching
         *     phi = 90 deg:  zero-crossing for debunching
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        ShortRF (
            amrex::ParticleReal V,
            amrex::ParticleReal freq,
            amrex::ParticleReal phase,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_V(V), m_freq(freq), m_phase(phase)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a shortrf functor, so that a variable of this type can be used like a
         *  shortrf function.
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

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            amrex::ParticleReal const k = (2.0_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);

            // access reference particle values (final, initial):
            amrex::ParticleReal const ptf_ref = refpart.pt;
            amrex::ParticleReal const pti_ref = ptf_ref + m_V*cos(phi);
            amrex::ParticleReal const bgf = sqrt(pow(ptf_ref, 2) - 1.0_prt);
            amrex::ParticleReal const bgi = sqrt(pow(pti_ref, 2) - 1.0_prt);

            // initial conversion from static to dynamic units:
            px = px*bgi;
            py = py*bgi;
            pt = pt*bgi;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum in dynamic units
            p.pos(RealAoS::x) = x;
            pxout = px;

            p.pos(RealAoS::y) = y;
            pyout = py;

            p.pos(RealAoS::t) = t;
            ptout = pt - m_V*cos(k*t + phi) + m_V*cos(phi);

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // final conversion from dynamic to static units:
            px = px/bgf;
            py = py/bgf;
            pt = pt/bgf;

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);

            // compute intial value of beta*gamma
            amrex::ParticleReal const bgi = sqrt(pow(pt, 2) - 1.0_prt);

            // advance pt
            refpart.pt = pt - m_V*cos(phi);

            // compute final value of beta*gamma
            amrex::ParticleReal const ptf = refpart.pt;
            amrex::ParticleReal const bgf = sqrt(pow(ptf, 2) - 1.0_prt);

            // advance position (x,y,z,t)
            refpart.x = x;
            refpart.y = y;
            refpart.z = z;
            refpart.t = t;

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

        }


    private:
        amrex::ParticleReal m_V; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_freq; //! RF frequency in Hz.
        amrex::ParticleReal m_phase; //! reference RF phase in degrees.
    };

} // namespace impactx

#endif // IMPACTX_SHORTRF_H
