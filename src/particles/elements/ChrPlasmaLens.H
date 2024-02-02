/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRPLASMALENS_H
#define IMPACTX_CHRPLASMALENS_H

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
    struct ChrPlasmaLens
    : public elements::BeamOptic<ChrPlasmaLens>,
      public elements::Thick,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "ChrPlasmaLens";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An active cylindrically symmetric plasma lens with chromatic focusing
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m.
         * @param g Azimuthal magnetic field gradient in T/m (g > 0)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param nslice number of slices used for the application of space charge
         */
        ChrPlasmaLens (
            amrex::ParticleReal ds,
            amrex::ParticleReal g,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int nslice = 1
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          m_g(g)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a chrquad functor, so that a variable of this type can be used like a chrquad function.
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

            // access reference particle values to find beta
            amrex::ParticleReal const bet = refpart.beta();

            // normalize focusing strength to magnetic rigidity
            amrex::ParticleReal g = m_g / refpart.rigidity_Tm();

            // compute particle momentum deviation delta + 1
            amrex::ParticleReal delta1;
            delta1 = sqrt(1_prt - 2_prt*pt/bet + pow(pt,2));
            amrex::ParticleReal const delta = delta1 - 1_prt;

            // compute phase advance per unit length in s (in rad/m)
            // chromatic dependence on delta is included
            amrex::ParticleReal const omega = sqrt(std::abs(g)/delta1);

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal const ptout = pt;

            // paceholder variables
            amrex::ParticleReal q1 = x;
            amrex::ParticleReal q2 = y;
            amrex::ParticleReal p1 = px;
            amrex::ParticleReal p2 = py;

            // advance transverse position and momentum (focusing quad)
            p.pos(RealAoS::x) = cos(omega*slice_ds)*x +
                                sin(omega*slice_ds)/(omega*delta1)*px;
            pxout = -omega*delta1*sin(omega*slice_ds)*x + cos(omega*slice_ds)*px;

            p.pos(RealAoS::y) = cos(omega*slice_ds)*y +
                                sin(omega*slice_ds)/(omega*delta1)*py;
            pyout = -omega*delta1*sin(omega*slice_ds)*y + cos(omega*slice_ds)*py;

            // advance longitudinal position and momentum

            // the corresponding symplectic update to t
            amrex::ParticleReal const term = pt + delta/bet;
            amrex::ParticleReal const t0 = t - term*slice_ds/delta1;

            amrex::ParticleReal const w = omega*delta1;
            amrex::ParticleReal const term1 = -(pow(p2,2)-pow(q2,2)*pow(w,2))*sin(2_prt*slice_ds*omega);
            amrex::ParticleReal const term2 = -(pow(p1,2)-pow(q1,2)*pow(w,2))*sin(2_prt*slice_ds*omega);
            amrex::ParticleReal const term3 = -2_prt*q2*p2*w*cos(2_prt*slice_ds*omega);
            amrex::ParticleReal const term4 = -2_prt*q1*p1*w*cos(2_prt*slice_ds*omega);
            amrex::ParticleReal const term5 = 2_prt*omega*(q1*p1*delta1 + q2*p2*delta1
                                        -(pow(p1,2)+pow(p2,2))*slice_ds - (pow(q1,2)+pow(q2,2))*pow(w,2)*slice_ds);
            p.pos(RealAoS::t) = t0 + (-1_prt+bet*pt)/(8_prt*bet*pow(delta1,3)*omega)
                                        *(term1+term2+term3+term4+term5);

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

        amrex::ParticleReal m_g; //! magnetic field gradient in T/m
    };

} // namespace impactx

#endif // IMPACTX_CHRPLASMALENS_H
