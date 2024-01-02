/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_MULTIPOLE_H
#define IMPACTX_MULTIPOLE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

namespace impactx
{
    struct Multipole
    : public elements::BeamOptic<Multipole>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "Multipole";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A general thin multipole element
         *
         * @param multipole index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
         * @param K_normal Integrated normal multipole coefficient (1/meter^m)
         * @param K_skew Integrated skew multipole coefficient (1/meter^m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        Multipole (
            int multipole,
            amrex::ParticleReal K_normal,
            amrex::ParticleReal K_skew,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_multipole(multipole), m_Kn(K_normal), m_Ks(K_skew)
        {
            // compute factorial of multipole index
            int const m = m_multipole - 1;
            m_mfactorial = 1;
            for( int n = 1; n < m + 1; n = n + 1 ) {
               m_mfactorial *= n;
            }
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
         * @param refpart reference particle (unused)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                PType& AMREX_RESTRICT p,
                amrex::ParticleReal & AMREX_RESTRICT px,
                amrex::ParticleReal & AMREX_RESTRICT py,
                amrex::ParticleReal & AMREX_RESTRICT pt,
                [[maybe_unused]] RefPart const & refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // a complex type with two amrex::ParticleReal
            using Complex = amrex::GpuComplex<amrex::ParticleReal>;

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // access reference particle values to find (beta*gamma)^2
            //amrex::ParticleReal const pt_ref = refpart.pt;
            //amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // assign complex position and complex multipole strength
            Complex const zeta(x, y);
            Complex const alpha(m_Kn, m_Ks);

            // compute complex momentum kick
            int const m = m_multipole - 1;
            Complex kick = amrex::pow(zeta, m);
            kick *= alpha;
            amrex::ParticleReal const dpx = -1.0_prt*kick.m_real/m_mfactorial;
            amrex::ParticleReal const dpy = kick.m_imag/m_mfactorial;

            // advance position and momentum
            p.pos(RealAoS::x) = x;
            pxout = px + dpx;

            p.pos(RealAoS::y) = y;
            pyout = py + dpy;

            p.pos(RealAoS::t) = t;
            ptout = pt;

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
        int m_multipole; //! multipole index
        int m_mfactorial; //! factorial of multipole index
        amrex::ParticleReal m_Kn; //! integrated normal multipole coefficient
        amrex::ParticleReal m_Ks; //! integrated skew multipole coefficient

    };

} // namespace impactx

#endif // IMPACTX_MULTIPOLE_H
