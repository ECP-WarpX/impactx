/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SBEND_H
#define IMPACTX_SBEND_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct Sbend
    {
        static constexpr auto name = "Sbend";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An ideal sector bend
         *
         * @param ds Segment length in m.
         * @param rc Radius of curvature in m.
         * @param nslice number of slices used for the application of space charge
         */
        Sbend( amrex::ParticleReal const ds, amrex::ParticleReal const rc,
               int const nslice)
        : m_ds(ds), m_rc(rc), m_nslice(nslice)
        {
        }

        /** This is a sbend functor, so that a variable of this type can be used like a sbend function.
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

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;
            amrex::ParticleReal const bet = sqrt(betgam2/(1.0_prt + betgam2));

            // calculate expensive terms once
            //   TODO: use sincos function once wrapped in AMReX
            amrex::ParticleReal const theta = slice_ds/m_rc;
            amrex::ParticleReal const sin_theta = sin(theta);
            amrex::ParticleReal const cos_theta = cos(theta);

            // advance position and momentum (sector bend)
            p.pos(0) = cos_theta*x + m_rc*sin_theta*px
                       - (m_rc/bet)*(1.0_prt - cos_theta)*pt;

            pxout = -sin_theta/m_rc*x + cos_theta*px - sin_theta/bet*pt;

            p.pos(1) = y + m_rc*theta*py;

            // pyout = py;

            p.pos(2) = sin_theta/bet*x + m_rc/bet*(1.0_prt - cos_theta)*px + t
                       + m_rc*(-theta+sin_theta/(bet*bet))*pt;

            // ptout = pt;

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
            amrex::ParticleReal const s = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // assign intermediate parameter
            amrex::ParticleReal const theta = slice_ds/m_rc;
            amrex::ParticleReal const B = sqrt(pow(pt,2)-1.0_prt)/m_rc;

            // calculate expensive terms once
            //   TODO: use sincos function once wrapped in AMReX
            amrex::ParticleReal const sin_theta = sin(theta);
            amrex::ParticleReal const cos_theta = cos(theta);

            // advance position and momentum (bend)
            refpart.px = px*cos_theta - pz*sin_theta;
            refpart.py = py;
            refpart.pz = pz*cos_theta + px*sin_theta;
            refpart.pt = pt;

            refpart.x = x + (refpart.pz - pz)/B;
            refpart.y = y + (theta/B)*py;
            refpart.z = z - (refpart.px - px)/B;
            refpart.t = t - (theta/B)*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;

        }

        /** Number of slices used for the application of space charge
         *
         * @return positive integer
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return m_nslice;
        }

        /** Return the segment length
         *
         * @return value in meters
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            return m_ds;
        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_rc; //! bend radius in m
        int m_nslice; //! number of slices used for the application of space charge
    };

} // namespace impactx

#endif // IMPACTX_SBEND_H
