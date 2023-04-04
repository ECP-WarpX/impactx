/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PROT_H
#define IMPACTX_PROT_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct PRot
    : public elements::BeamOptic<PRot>,
      public elements::Thin
    {
        static constexpr auto name = "PRot";
        using PType = ImpactXParticleContainer::ParticleType;

        constexpr amrex::ParticleReal degree2rad = ablastr::constant::math::pi / 180.0;

        /** An exact pole face rotation in the x-z plane, from a frame
         *  in which the reference orbit has angle phi_in with the z-axis,
         *  to a frame in which the reference orbit has angle phi_out
         *  with the z-axis.  The net rotation angle is phi_out - phi_in.
         *
         * @param phi_in Initial angle of reference trajectory w/r/t z (degrees)
         * @param phi_out Final angle of reference trajectory w/r/t/ z (degrees)
         */
        PRot( amrex::ParticleReal const phi_in, amrex::ParticleReal const phi_out )
        : m_phi_in(phi_in * degree2rad), m_phi_out(phi_out * degree2rad)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a prot functor, so that a variable of this type can be used like a
         *  prot function.
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

            // access reference particle values to find beta:
            amrex::ParticleReal const beta = refpart.beta();

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // store rotation angle and initial, final values of pz
            amrex::ParticleReal theta = m_phi_out - m_phi_in;
            amrex::ParticleReal pz = sqrt(1.0_prt - 2.0_prt*pt/beta
               + pow(pt,2) - pow(py,2) - pow(px + sin(m_phi_in),2));
            amrex::ParticleReal pzf = pz*cos(theta) - (px +
                 sin(m_phi_in))*sin(theta);

            // advance position and momentum
            p.pos(0) = x*pz/pzf;
            pxout = px*cos(theta) + (pz - cos(m_phi_in))*sin(theta);

            p.pos(1) = y + py*x*sin(theta)/pzf;
            pyout = py;

            p.pos(2) = t - (pt - 1.0_prt/beta)*x*sin(theta)/pzf;
            ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

        }

        /** This pushes the reference particle. */
        using Thin::operator();

    private:
        amrex::ParticleReal m_phi_in; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_phi_out; //! RF wavenumber in 1/m.
    };

} // namespace impactx

#endif // IMPACTX_PROT_H
