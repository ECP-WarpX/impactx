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
#include "mixin/regular.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

namespace impactx
{
    struct NonlinearLens
    : public elements::Regular<NonlinearLens>
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
         */
        NonlinearLens( amrex::ParticleReal const knll,
                       amrex::ParticleReal const cnll )
        : m_knll(knll), m_cnll(cnll)
        {
        }

        /** Push all particles */
        using Regular::operator();

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

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // access reference particle values to find (beta*gamma)^2
            //amrex::ParticleReal const pt_ref = refpart.pt;
            //amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // assign complex position zeta = x + iy
            Complex zeta(x, y);
            Complex re1(1.0_prt, 0.0_prt);
            Complex im1(0.0_prt, 1.0_prt);

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
            amrex::ParticleReal kick = -m_knll/m_cnll;
            amrex::ParticleReal dpx = kick*dF.m_real;
            amrex::ParticleReal dpy = -kick*dF.m_imag;

            // advance position and momentum
            p.pos(0) = x;
            pxout = px + dpx;

            p.pos(1) = y;
            pyout = py + dpy;

            p.pos(2) = t;
            ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() ([[maybe_unused]] RefPart & AMREX_RESTRICT refpart) const {

            // nothing to do: this is a zero-length element
        }

        /** Number of slices used for the application of space charge
         *
         * @return one, because this is a zero-length element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return 1;
        }

        /** Return the segment length
         *
         * @return zero, because this is a zero-length element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            using namespace amrex::literals;
            return 0.0_prt;
        }

    private:
        amrex::ParticleReal m_knll; //! integrated strength of the nonlinear lens (m)
        amrex::ParticleReal m_cnll; //! distance of singularities from the origin (m)
    };

} // namespace impactx

#endif // IMPACTX_NONLINEARLENS_H
