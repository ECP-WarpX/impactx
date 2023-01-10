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
#include "mixin/regular.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct ShortRF
    : public elements::Regular<ShortRF>
    {
        static constexpr auto name = "ShortRF";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A short RF cavity element at zero crossing for bunching
         *
         * @param V Normalized RF voltage drop V = Emax*L/(c*Brho)
         * @param k Wavenumber of RF in 1/m
         */
        ShortRF( amrex::ParticleReal const V, amrex::ParticleReal const k )
        : m_V(V), m_k(k)
        {
        }

        /** Push all particles */
        using Regular::operator();

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

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // access reference particle values to find (beta*gamma)^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum
            p.pos(0) = x;
            pxout = px + m_k*m_V/(2.0_prt*betgam2)*x;

            p.pos(1) = y;
            pyout = py + m_k*m_V/(2.0_prt*betgam2)*y;

            p.pos(2) = t;
            ptout = pt - m_k*m_V*t;

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
        amrex::ParticleReal m_V; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_k; //! RF wavenumber in 1/m.
    };

} // namespace impactx

#endif // IMPACTX_SHORTRF_H
