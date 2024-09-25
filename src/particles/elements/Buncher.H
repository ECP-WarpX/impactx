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
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct Buncher
    : public elements::Named,
      public elements::BeamOptic<Buncher>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto type = "Buncher";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A short RF cavity element at zero crossing for bunching
         *
         * @param V Normalized RF voltage drop V = Emax*L/(c*Brho)
         * @param k Wavenumber of RF in 1/m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        Buncher (
            amrex::ParticleReal V,
            amrex::ParticleReal k,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(name),
          Alignment(dx, dy, rotation_degree),
          m_V(V), m_k(k)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a buncher functor, so that a variable of this type can be used like a
         *  buncher function.
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
                RefPart const & refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // access reference particle values to find (beta*gamma)^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum
            pxout = px + m_k*m_V/(2.0_prt*betgam2)*x;
            pyout = py + m_k*m_V/(2.0_prt*betgam2)*y;
            ptout = pt - m_k*m_V*t;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        amrex::ParticleReal m_V; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_k; //! RF wavenumber in 1/m.
    };

} // namespace impactx

#endif // IMPACTX_BUNCHER_H
