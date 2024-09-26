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
#include "mixin/alignment.H"
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
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto type = "ThinDipole";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal degree2rad = ablastr::constant::math::pi / 180.0;

        /** A general thin-kick dipole element with chromatic effects
         *
         * This model is equivalent to the model described in Section 3.1 of
         * G. Ripken, F. Schmidt, "A Symplectic Six-Dimensional Thin-Lens Formalism
         * for Tracking," CERN/SL/95-12 (AP), 1995.  It is intended to replicate
         * the thin-lens dipole model in MAD-X.
         *
         * @param theta - the total bending angle (degrees)
         * @param rc - the curvature radius (m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        ThinDipole (
            amrex::ParticleReal theta,
            amrex::ParticleReal rc,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_theta(theta * degree2rad),
          m_rc(rc)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a multipole functor, so that a variable of this type can be used like a
         *  multipole function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] uint64_t & AMREX_RESTRICT idcpu,
            RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // access position data
            amrex::ParticleReal const xout = x;
            amrex::ParticleReal const yout = y;
            amrex::ParticleReal const tout = t;

            // access reference particle to find relativistic beta
            amrex::ParticleReal const beta_ref = refpart.beta();

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // compute the function expressing dp/p in terms of pt (labeled f in Ripken etc.)
            amrex::ParticleReal f = -1.0_prt + std::sqrt(1.0_prt - 2.0_prt*pt/beta_ref + std::pow(pt,2));
            amrex::ParticleReal fprime = (1.0_prt - beta_ref*pt)/(beta_ref*(1.0_prt + f));

            // compute the effective (equivalent) arc length and curvature
            amrex::ParticleReal ds = m_theta*m_rc;
            amrex::ParticleReal kx = 1.0_prt/m_rc;

            // advance position and momentum
            x = xout;
            pxout = px - std::pow(kx,2) * ds * xout + kx * ds * f; //eq (3.2b)

            y = yout;
            pyout = py;

            t = tout + kx * xout * ds * fprime; //eq (3.2e)
            ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            // assign input reference particle values
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const pz = refpart.pz;

            // calculate expensive terms once
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(m_theta);

            // advance position and momentum (thin dipole)
            refpart.px = px*cos_theta - pz*sin_theta;
            refpart.pz = pz*cos_theta + px*sin_theta;
        }

        amrex::ParticleReal m_theta; //! dipole bending angle (rad)
        amrex::ParticleReal m_rc; //! curvature radius (m)

    };

} // namespace impactx

#endif // IMPACTX_THINDIPOLE_H
