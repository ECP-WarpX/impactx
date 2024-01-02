/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_QUAD_H
#define IMPACTX_QUAD_H

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
    struct Quad
    : public elements::BeamOptic<Quad>,
      public elements::Thick,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "Quad";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A Quadrupole magnet
         *
         * @param ds Segment length in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param nslice number of slices used for the application of space charge
         */
        Quad (
            amrex::ParticleReal ds,
            amrex::ParticleReal k,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int nslice = 1
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          m_k(k)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a quad functor, so that a variable of this type can be used like a quad function.
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
            RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = sqrt(std::abs(m_k));

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal const ptout = pt;

            if(m_k > 0.0) {
               // advance position and momentum (focusing quad)
               p.pos(RealAoS::x) = cos(omega*slice_ds)*x + sin(omega*slice_ds)/omega*px;
               pxout = -omega*sin(omega*slice_ds)*x + cos(omega*slice_ds)*px;

               p.pos(RealAoS::y) = cosh(omega*slice_ds)*y + sinh(omega*slice_ds)/omega*py;
               pyout = omega*sinh(omega*slice_ds)*y + cosh(omega*slice_ds)*py;

               p.pos(RealAoS::t) = t + (slice_ds/betgam2)*pt;
               // ptout = pt;
            } else {
               // advance position and momentum (defocusing quad)
               p.pos(RealAoS::x) = cosh(omega*slice_ds)*x + sinh(omega*slice_ds)/omega*px;
               pxout = omega*sinh(omega*slice_ds)*x + cosh(omega*slice_ds)*px;

               p.pos(RealAoS::y) = cos(omega*slice_ds)*y + sin(omega*slice_ds)/omega*py;
               pyout = -omega*sin(omega*slice_ds)*y + cos(omega*slice_ds)*py;

               p.pos(RealAoS::t) = t + (slice_ds/betgam2)*pt;
               // ptout = pt;
            }

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
        amrex::ParticleReal m_k; //! quadrupole strength in 1/m
    };

} // namespace impactx

#endif // IMPACTX_QUAD_H
