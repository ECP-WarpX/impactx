/* Copyright 2022 The Regents of the University of California, through Lawrence
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

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

namespace impactx
{
    struct Multipole
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A general thin multipole element
         *
         * @param multipole index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
         * @param K_normal Integrated normal multipole coefficient (1/meter^m)
         * @param K_skew Integrated skew multipole coefficient (1/meter^m)
         */
        Multipole( int const multipole,
                   amrex::ParticleReal const K_normal,
                   amrex::ParticleReal const K_skew )
        : m_multipole(multipole), m_Kn(K_normal), m_Ks(K_skew)
        {
            // compute factorial of multipole index
            int const m = m_multipole - 1;
            m_mfactorial = 1;
            for( int n = 1; n < m + 1; n = n + 1 ) {
               m_mfactorial *= n;
            }
        }

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
        void operator() (
                [[maybe_unused]] RefPart & AMREX_RESTRICT refpart) const {

            // nothing to do: this is a zero-length element
        }

        /** Number of slices used for the application of space charge
         *
         * @return zero, because this is a zero-length element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return 0;
        }

    private:
        int m_multipole; //! multipole index
        int m_mfactorial; //! factorial of multipole index
        amrex::ParticleReal m_Kn; //! integrated normal multipole coefficient
        amrex::ParticleReal m_Ks; //! integrated skew multipole coefficient

    };

} // namespace impactx

#endif // IMPACTX_MULTIPOLE_H
