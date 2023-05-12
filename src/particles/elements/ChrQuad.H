/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRQUAD_H
#define IMPACTX_CHRQUAD_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct ChrQuad
    : public elements::BeamOptic<ChrQuad>,
      public elements::Thick,
      public elements::NoFinalize
    {
        static constexpr auto name = "ChrQuad";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A Quadrupole magnet with chromatic focusing
         *
         * @param ds Segment length in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param nslice number of slices used for the application of space charge
         */
        ChrQuad( amrex::ParticleReal const ds, amrex::ParticleReal const k,
              int const nslice )
        : Thick(ds, nslice), m_k(k)
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
                RefPart const & refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta, gamma
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const gam = refpart.gamma();

            // compute particle momentum deviation delta + 1
            amrex::ParticleReal delta1;
            delta1 = sqrt(1_prt - 2_prt*pt/bet + pow(pt,2));
            amrex::ParticleReal delta = delta1 - 1_prt;

            // compute phase advance per unit length in s (in rad/m)
            // chromatic dependence on delta is included
            amrex::ParticleReal const omega = sqrt(std::abs(m_k)/delta1);

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            if(m_k > 0.0) {
               // advance transverse position and momentum (focusing quad)
               p.pos(RealAoS::x) = cos(omega*slice_ds)*x +
                                   sin(omega*slice_ds)/(omega*delta1)*px;
               pxout = -omega*delta1*sin(omega*slice_ds)*x + cos(omega*slice_ds)*px;

               p.pos(RealAoS::y) = cosh(omega*slice_ds)*y +
                                   sinh(omega*slice_ds)/(omega*delta1)*py;
               pyout = omega*delta1*sinh(omega*slice_ds)*y + cosh(omega*slice_ds)*py;

            } else {
               // advance transverse position and momentum (defocusing quad)
               p.pos(RealAoS::x) = cosh(omega*slice_ds)*x +
                                   sinh(omega*slice_ds)/(omega*delta1)*px;
               pxout = omega*delta1*sinh(omega*slice_ds)*x + cosh(omega*slice_ds)*px;

               p.pos(RealAoS::y) = cos(omega*slice_ds)*y +
                                   sin(omega*slice_ds)/(omega*delta1)*py;
               pyout = -omega*delta1*sin(omega*slice_ds)*y + cos(omega*slice_ds)*py;

            }

            // advance longitudinal position and momentum

            // the corresponding symplectic update to t
            amrex::ParticleReal term = pt + delta/bet;
        amrex::ParticleReal t0 = t - term*slice_ds/delta1;

        amrex::ParticleReal w = omega*delta1;
            amrex::ParticleReal term1 = -(pow(py,2)+pow(y,2)*pow(w,2))*sinh(2_prt*slice_ds*omega);
        amrex::ParticleReal term2 = -(pow(px,2)-pow(x,2)*pow(w,2))*sin(2_prt*slice_ds*omega);
        amrex::ParticleReal term3 = -2_prt*y*py*w*cosh(2_prt*slice_ds*omega);
        amrex::ParticleReal term4 = -2_prt*x*px*w*cos(2_prt*slice_ds*omega);
            amrex::ParticleReal term5 = 2_prt*omega*(x*px*delta1 + y*py*delta1
                -(pow(px,2)+pow(py,2))*slice_ds - (pow(x,2)-pow(y,2))*pow(w,2)*slice_ds);
            p.pos(RealAoS::t) = t0 + (-1_prt+bet*pt)/(8_prt*bet*pow(delta1,3)*omega)
                                        *(term1+term2+term3+term4+term5);

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

#endif // IMPACTX_CHRQUAD_H
