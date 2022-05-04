/* Copyright 2022 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_NONLINEARLENS_H
#define IMPACTX_NONLINEARLENS_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

namespace impactx
{
    struct NonlinearLens
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A thin lens associated with a single short segment of the
         *  nonlinear magnetic insert described by V. Danilov and
         *  S. Nagaitsev, PRSTAB 13, 084002 (2010), Sect. V.A.  This
         *  element appears in MAD-X as type NLLENS.
         *
         * @param knll - integrated strength of the nonlinear lens (m)
         * @param cnll - distance of singularities from the origin (m)
         *
         */
        NonlinearLens( amrex::ParticleReal const knll,
                   amrex::ParticleReal const cnll )
        : m_knll(knll), m_cnll(cnll)
        {
        }

        /** This is a nonlinear lens functor, so that a variable of this type can be used like a
         *  nonlinear lens function.
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
                [[maybe_unused]] RefPart const refpart) const {

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
            amrex::GpuComplex zeta(x,y);
            amrex::GpuComplex re1(1.0,0.0);
            amrex::GpuComplex im1(0.0,1.0);

            // compute croot = sqrt(1-zeta**2)
            amrex::GpuComplex croot = amrex::pow(zeta,2);
            croot = re1 - croot;
            croot = amrex::sqrt(croot);

            // compute carcsin = arcsin(zeta)
            amrex::GpuComplex carcsin = im1*zeta + croot;
            carcsin = -im1*amrex::log(carcsin);

            // compute complex function F'(zeta)
            amrex::GpuComplex dF = zeta/amrex::pow(croot,2);
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

    private:
        amrex::ParticleReal m_knll; //! integrated strength of the nonlinear lens (m)
        amrex::ParticleReal m_cnll; //! distance of singularities from the origin (m)

    };

} // namespace impactx

#endif // IMPACTX_NONLINEARLENS_H
