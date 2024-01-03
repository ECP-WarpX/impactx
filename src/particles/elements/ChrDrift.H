/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRDRIFT_H
#define IMPACTX_CHRDRIFT_H

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
    struct ChrDrift
    : public elements::BeamOptic<ChrDrift>,
      public elements::Thick,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "ChrDrift";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A drift with chromatic effects included
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param nslice number of slices used for the application of space charge
         */
        ChrDrift (
            amrex::ParticleReal ds,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int nslice = 1
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a chrdrift functor, so that a variable of this type can be used like a chrdrift function.
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

            // initialize output values of momenta
            amrex::ParticleReal const pxout = px;
            amrex::ParticleReal const pyout = py;
            amrex::ParticleReal const ptout = pt;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta, gamma
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const gam = refpart.gamma();

            // compute particle momentum deviation delta + 1
            amrex::ParticleReal delta1;
            delta1 = sqrt(1_prt - 2_prt*pt/bet + pow(pt,2));

            // advance transverse position and momentum (drift)
            p.pos(RealAoS::x) = x + slice_ds * px / delta1;
            // pxout = px;
            p.pos(RealAoS::y) = y + slice_ds * py / delta1;
            // pyout = py;

            // the corresponding symplectic update to t
            amrex::ParticleReal term = 2_prt*pow(pt,2)+pow(px,2)+pow(py,2);
            term = 2_prt - 4_prt*bet*pt + pow(bet,2)*term;
            term = -2_prt + pow(gam,2)*term;
            term = (-1_prt+bet*pt)*term;
            term = term/(2_prt*pow(bet,3)*pow(gam,2));
            p.pos(RealAoS::t) = t - slice_ds*(1_prt/bet + term/pow(delta1,3));
            // ptout = pt;

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

            // advance position and momentum (drift)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;

        }
    };

} // namespace impactx

#endif // IMPACTX_CHRDRIFT_H
