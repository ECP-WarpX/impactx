/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CFBEND_H
#define IMPACTX_CFBEND_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct CFbend
    : public elements::BeamOptic<CFbend>,
      public elements::Thick,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "CFbend";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An combined-function bend, consisting of an ideal sector bend with
         *  an upright quadrupole focusing component.
         *
         * @param ds Segment length in m.
         * @param rc Radius of curvature in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param nslice number of slices used for the application of space charge
         */
        CFbend (
            amrex::ParticleReal ds,
            amrex::ParticleReal rc,
            amrex::ParticleReal k,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int nslice = 1
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          m_rc(rc), m_k(k)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a cfbend functor, so that a variable of this type can be used like a cfbend function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] uint64_t & AMREX_RESTRICT idcpu,
            RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // initialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;

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

            // update horizontal and longitudinal phase space variables
            amrex::ParticleReal const gx = m_k + pow(m_rc,-2);
            amrex::ParticleReal const omegax = sqrt(std::abs(gx));

            if(gx > 0.0) {
                // calculate expensive terms once
                auto const [sinx, cosx] = amrex::Math::sincos(omegax * slice_ds);
                amrex::ParticleReal const r56 = slice_ds/betgam2
                    + (sinx - omegax*slice_ds)/(gx*omegax*pow(bet,2)*pow(m_rc,2));

                // advance position and momentum (focusing)
                x = cosx*xout + sinx/omegax*px - (1.0_prt - cosx)/(gx*bet*m_rc)*pt;
                pxout = -omegax*sinx*xout + cosx*px - sinx/(omegax*bet*m_rc)*pt;

                y = sinx/(omegax*bet*m_rc)*xout + (1.0_prt - cosx)/(gx*bet*m_rc)*px
                    + tout + r56*pt;
                ptout = pt;
            } else {
                // calculate expensive terms once
                amrex::ParticleReal const sinhx = sinh(omegax * slice_ds);
                amrex::ParticleReal const coshx = cosh(omegax * slice_ds);
                amrex::ParticleReal const r56 = slice_ds/betgam2
                    + (sinhx - omegax*slice_ds)/(gx*omegax*pow(bet,2)*pow(m_rc,2));

                // advance position and momentum (defocusing)
                x = coshx*xout + sinhx/omegax*px - (1.0_prt - coshx)/(gx*bet*m_rc)*pt;
                pxout = omegax*sinhx*xout + coshx*px - sinhx/(omegax*bet*m_rc)*pt;

                t = sinhx/(omegax*bet*m_rc)*xout + (1.0_prt - coshx)/(gx*bet*m_rc)*px
                    + tout + r56*pt;
                ptout = pt;
            }

            // update vertical phase space variables
            amrex::ParticleReal const gy = -m_k;
            amrex::ParticleReal const omegay = sqrt(std::abs(gy));

            if(gy > 0.0) {
                // calculate expensive terms once
                auto const [siny, cosy] = amrex::Math::sincos(omegay * slice_ds);

                // advance position and momentum (focusing)
                y = cosy*yout + siny/omegay*py;
                pyout = -omegay*siny*yout + cosy*py;

            } else {
                // calculate expensive terms once
                amrex::ParticleReal const sinhy = sinh(omegay * slice_ds);
                amrex::ParticleReal const coshy = cosh(omegay * slice_ds);

                // advance position and momentum (defocusing)
                y = coshy*yout + sinhy/omegay*py;
                pyout = omegay*sinhy*yout + coshy*py;
            }

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
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
            amrex::ParticleReal const theta = slice_ds/m_rc;
            amrex::ParticleReal const B = sqrt(pow(pt,2)-1.0_prt)/m_rc;

            // calculate expensive terms once
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

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

        amrex::ParticleReal m_rc; //! bend radius in m
        amrex::ParticleReal m_k;  //! quadrupole strength in m^(-2)
    };

} // namespace impactx

#endif // IMPACTX_CFBEND_H
