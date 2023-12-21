/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRACC_H
#define IMPACTX_CHRACC_H

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
    struct ChrAcc
    : public elements::BeamOptic<ChrAcc>,
      public elements::Thick,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "ChrAcc";
        using PType = ImpactXParticleContainer::ParticleType;

        /** Acceleration in a uniform field Ez, with a uniform solenoidal field Bz.
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m.
         * @param ez electric field strength in m^(-1)
         *           = (charge * electric field Ez in V/m) / (m*c^2)
         * @param bz magnetic field strength in m^(-1)
         *           = (charge * magnetic field Bz in T) / (m*c)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error around the s axis (in the x-y plane) [degrees]
         * @param nslice number of slices used for the application of space charge
         */
        ChrAcc (
            amrex::ParticleReal ds,
            amrex::ParticleReal ez,
            amrex::ParticleReal bz,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int nslice = 1
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          m_ez(ez), m_bz(bz)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a chracc functor, so that a variable of this type can be used like a chracc function.
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

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values (final, initial):
            amrex::ParticleReal const ptf_ref = refpart.pt;
            amrex::ParticleReal const pti_ref = ptf_ref + m_ez*slice_ds;
            amrex::ParticleReal const bgf = sqrt(pow(ptf_ref, 2) - 1.0_prt);
            amrex::ParticleReal const bgi = sqrt(pow(pti_ref, 2) - 1.0_prt);

            // initial conversion from static to dynamic units:
            px = px*bgi;
            py = py*bgi;
            pt = pt*bgi;

            // compute intermediate quantities related to acceleration
            amrex::ParticleReal const pti_tot = pti_ref + pt;
            amrex::ParticleReal const ptf_tot = ptf_ref + pt;
            amrex::ParticleReal const pzi_tot = sqrt(pow(pti_tot,2)-1_prt);
        amrex::ParticleReal const pzf_tot = sqrt(pow(ptf_tot,2)-1_prt);
            amrex::ParticleReal const pzi_ref = sqrt(pow(pti_ref,2)-1_prt);
        amrex::ParticleReal const pzf_ref = sqrt(pow(ptf_ref,2)-1_prt);

        amrex::ParticleReal const numer = -ptf_tot + pzf_tot;
            amrex::ParticleReal const denom = -pti_tot + pzi_tot;

            // compute focusing constant (1/m) and rotation angle (in rad)
            amrex::ParticleReal const alpha = m_bz/2.0_prt;
            amrex::ParticleReal const theta = alpha/m_ez*log(numer/denom);

            // intialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance positions and momenta using map for focusing
            xout = cos(theta)*x + sin(theta)/alpha*px;
            pxout = -alpha*sin(theta)*x + cos(theta)*px;

            yout = cos(theta)*y + sin(theta)/alpha*py;
            pyout = -alpha*sin(theta)*y + cos(theta)*py;

            // the correct symplectic update for t
            tout = t + (pzf_tot - pzf_ref - pzi_tot + pzi_ref)/m_ez;
        tout = tout + (1_prt/pzi_tot - 1_prt/pzf_tot)*(pow(py-alpha*x,2)+pow(px+alpha*y,2))/(2_prt*m_ez);
        ptout = pt;

            // assign intermediate momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // advance positions and momenta using map for rotation
            p.pos(RealAoS::x) = cos(theta)*xout + sin(theta)*yout;
            pxout = cos(theta)*px + sin(theta)*py;

            p.pos(RealAoS::y) = -sin(theta)*xout + cos(theta)*yout;
            pyout = -sin(theta)*px + cos(theta)*py;

            p.pos(RealAoS::t) = tout;
            ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // final conversion from dynamic to static units:
            px = px/bgf;
            py = py/bgf;
            pt = pt/bgf;

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

            // compute intial value of beta*gamma
            amrex::ParticleReal const bgi = sqrt(pow(pt, 2) - 1.0_prt);

            // advance pt (uniform acceleration)
            refpart.pt = pt - m_ez*slice_ds;

            // compute final value of beta*gamma
            amrex::ParticleReal const ptf = refpart.pt;
            amrex::ParticleReal const bgf = sqrt(pow(ptf, 2) - 1.0_prt);

            // update t
            refpart.t = t + (bgf - bgi)/m_ez;

            // advance position (x,y,z)
            refpart.x = x + slice_ds*px/bgi;
            refpart.y = y + slice_ds*py/bgi;
            refpart.z = z + slice_ds*pz/bgi;

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

    private:
        amrex::ParticleReal m_ez; //! electric field strength in 1/m
        amrex::ParticleReal m_bz; //! magnetic field strength in 1/m
    };

} // namespace impactx

#endif // IMPACTX_CHRACC_H
