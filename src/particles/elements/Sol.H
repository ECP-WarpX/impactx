/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOL_H
#define IMPACTX_SOL_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct Sol
    : public elements::BeamOptic<Sol>,
      public elements::Thick,
      public elements::NoFinalize
    {
        static constexpr auto name = "Sol";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An ideal hard-edge Solenoid magnet
         *
         * @param ds Segment length in m.
         * @param ks Solenoid strength in m^(-1) (MADX convention)
         *           = (magnetic field Bz in T) / (rigidity in T-m)
         * @param nslice number of slices used for the application of space charge
         */
        Sol (
            amrex::ParticleReal ds,
            amrex::ParticleReal ks,
            int nslice
        )
        : Thick(ds, nslice), m_ks(ks)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a sol functor, so that a variable of this type can be used like a sol function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            PType & AMREX_RESTRICT p,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // compute phase advance per unit length (in rad/m) and
            // rotation angle (in rad)
            amrex::ParticleReal const alpha = m_ks/2.0_prt;
            amrex::ParticleReal const theta = alpha*slice_ds;

            // intialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance positions and momenta using map for focusing
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);
            xout = cos_theta*x + sin_theta/alpha*px;
            pxout = -alpha*sin_theta*x + cos_theta*px;

            yout = cos_theta*y + sin_theta/alpha*py;
            pyout = -alpha*sin_theta*y + cos_theta*py;

            tout = t + (slice_ds/betgam2)*pt;
            ptout = pt;

            // assign intermediate momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // advance positions and momenta using map for rotation
            p.pos(RealAoS::x) = cos_theta*xout + sin_theta*yout;
            pxout = cos_theta*px + sin_theta*py;

            p.pos(RealAoS::y) = -sin_theta*xout + cos_theta*yout;
            pyout = -sin_theta*px + cos_theta*py;

            p.pos(RealAoS::t) = tout;
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
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
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
            amrex::ParticleReal const step = slice_ds / sqrt(pow(pt,2)-1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

    private:
        amrex::ParticleReal m_ks; //! solenoid strength in 1/m
    };

} // namespace impactx

#endif // IMPACTX_SOL_H
