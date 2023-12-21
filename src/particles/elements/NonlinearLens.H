/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_NONLINEARLENS_H
#define IMPACTX_NONLINEARLENS_H

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
    struct NonlinearLens
    : public elements::BeamOptic<NonlinearLens>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "NonlinearLens";
        using PType = ImpactXParticleContainer::ParticleType;

        /** Single short segment of the nonlinear magnetic insert element
         *
         *  A thin lens associated with a single short segment of the
         *  nonlinear magnetic insert described by V. Danilov and
         *  S. Nagaitsev, PRSTAB 13, 084002 (2010), Sect. V.A.  This
         *  element appears in MAD-X as type NLLENS.
         *
         * @param knll integrated strength of the nonlinear lens (m)
         * @param cnll distance of singularities from the origin (m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error around the s axis (in the x-y plane) [degrees]
         */
        NonlinearLens (
            amrex::ParticleReal knll,
            amrex::ParticleReal cnll,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_knll(knll), m_cnll(cnll)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a nonlinear lens functor, so that a variable of this type
         *  can be used like a nonlinear lens function.
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

            // a complex type with two amrex::ParticleReal
            using Complex = amrex::GpuComplex<amrex::ParticleReal>;

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

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

            // assign complex position zeta = (x + iy)/cnll
            Complex zeta(x, y);
            zeta = zeta/m_cnll;
            Complex const re1(1.0_prt, 0.0_prt);
            Complex const im1(0.0_prt, 1.0_prt);

            // compute croot = sqrt(1-zeta**2)
            Complex croot = amrex::pow(zeta, 2);
            croot = re1 - croot;
            croot = amrex::sqrt(croot);

            // compute carcsin = arcsin(zeta)
            Complex carcsin = im1*zeta + croot;
            carcsin = -im1*amrex::log(carcsin);

            // compute complex function F'(zeta)
            Complex dF = zeta/amrex::pow(croot, 2);
            dF = dF + carcsin/amrex::pow(croot,3);

            // compute momentum kick
            amrex::ParticleReal const kick = -m_knll/m_cnll;
            amrex::ParticleReal const dpx = kick*dF.m_real;
            amrex::ParticleReal const dpy = -kick*dF.m_imag;

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
        amrex::ParticleReal m_knll; //! integrated strength of the nonlinear lens (m)
        amrex::ParticleReal m_cnll; //! distance of singularities from the origin (m)
    };

} // namespace impactx

#endif // IMPACTX_NONLINEARLENS_H
