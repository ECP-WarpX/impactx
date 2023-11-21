/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_THINDIPOLE_H
#define IMPACTX_THINDIPOLE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

namespace impactx
{
    struct ThinDipole
    : public elements::BeamOptic<ThinDipole>,
      public elements::Thin,
      public elements::NoFinalize
    {
        static constexpr auto name = "ThinDipole";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A general thin-kick dipole element with chromatic effects
         *
         * This model is equivalent to the model described in Section 3.1 of
         * G. Ripken, F. Schmidt, "A Symplectic Six-Dimensional Thin-Lens Formalism
         * for Tracking," CERN/SL/95-12 (AP), 1995.  It is intended to replicate
         * the thin-lens dipole model in MAD-X.
         *
         * @param theta - the total bending angle (degrees)
         * @param rc - the curvature radius (m)
         */
        ThinDipole( amrex::ParticleReal const theta,
                    amrex::ParticleReal const rc)
        : m_theta(theta * degree2rad), m_rc(rc)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a multipole functor, so that a variable of this type can be used like a
         *  multipole function.
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
                [[maybe_unused]] RefPart const & refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // access reference particle to find relativistic beta
            amrex::ParticleReal const beta_ref = refpart.beta();

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // compute the function expressing dp/p in terms of pt (labeled f in Ripken etc.)
            amrex::ParticleReal f = -1.0_prt + sqrt(1.0_prt - 2.0_prt*pt/beta_ref + pow(pt,2));
            amrex::ParticleReal fprime = (1.0_prt + beta_ref*pt)/(beta_ref*(1.0_prt + f);

            // compute the effective (equivalent) arc length and curvature
            amrex::ParticleReal ds = m_theta*m_rc;
            amrex::ParticleReal kx = 1.0_prt/rc;

            // advance position and momentum
            p.pos(RealAoS::x) = x;
            pxout = px - pow(kx,2)*ds*x + kx*ds*f; //eq (3.2b)

            p.pos(RealAoS::y) = y;
            pyout = py;

            p.pos(RealAoS::t) = t - kx*x*ds*fprime; //eq (3.2e)
            ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

        }

        /** This pushes the reference particle. */
        using Thin::operator();

    private:
        amrex::ParticleReal m_theta; //! dipole bending angle (rad)
        amrex::ParticleReal m_rc; //! curvature radius (m)

    };

} // namespace impactx

#endif // IMPACTX_THINDIPOLE_H
