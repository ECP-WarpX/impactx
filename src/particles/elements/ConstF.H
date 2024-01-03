/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CONSTF_H
#define IMPACTX_CONSTF_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct ConstF
    : public elements::BeamOptic<ConstF>,
      public elements::Thick,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "ConstF";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A linear Constant Focusing element
         *
         * @param ds Segment length in m.
         * @param kx Focusing strength for x in 1/m.
         * @param ky Focusing strength for y in 1/m.
         * @param kt Focusing strength for t in 1/m.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param nslice number of slices used for the application of space charge
         */
        ConstF (
            amrex::ParticleReal ds,
            amrex::ParticleReal kx,
            amrex::ParticleReal ky,
            amrex::ParticleReal kt,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int nslice = 1
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          m_kx(kx), m_ky(ky), m_kt(kt)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This pushes a single particle, relative to the reference particle
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

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // advance position and momentum
            p.pos(RealAoS::x) = cos(m_kx*slice_ds)*x + sin(m_kx*slice_ds)/m_kx*px;
            pxout = -m_kx*sin(m_kx*slice_ds)*x + cos(m_kx*slice_ds)*px;

            p.pos(RealAoS::y) = cos(m_ky*slice_ds)*y + sin(m_ky*slice_ds)/m_ky*py;
            pyout = -m_ky*sin(m_ky*slice_ds)*y + cos(m_ky*slice_ds)*py;

            p.pos(RealAoS::t) = cos(m_kt*slice_ds)*t + sin(m_kt*slice_ds)/(betgam2*m_kt)*pt;
            ptout = -(m_kt*betgam2)*sin(m_kt*slice_ds)*t + cos(m_kt*slice_ds)*pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
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
            amrex::ParticleReal const step = slice_ds / sqrt(pow(pt, 2)-1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;

        }

        amrex::ParticleReal m_kx; //! focusing x strength in 1/m
        amrex::ParticleReal m_ky; //! focusing y strength in 1/m
        amrex::ParticleReal m_kt; //! focusing t strength in 1/m
    };

} // namespace impactx

#endif // IMPACTX_CONSTF_H
