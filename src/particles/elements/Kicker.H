/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_KICKER_H
#define IMPACTX_KICKER_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

namespace impactx
{
    struct Kicker
    : public elements::BeamOptic<Kicker>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "Kicker";
        using PType = ImpactXParticleContainer::ParticleType;

        enum UnitSystem
        {
            dimensionless = 0, //! in units of the magnetic rigidity of the reference particle
            Tm = 1  //! T-m
        };

        /** A thin transverse kicker element, inducing a transverse
         *  momentum displacement from the reference orbit.
         *
         * @param xkick Strength of horizontal kick
         * @param ykick Strength of vertical kick
         * @param unit units of xkick and ykick
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        Kicker (
            amrex::ParticleReal xkick,
            amrex::ParticleReal ykick,
            UnitSystem unit,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_xkick(xkick), m_ykick(ykick), m_unit(unit)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a transverse kicker functor, so that a variable of this type can be
         * used like a function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param refpart reference particle (unused)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            PType& AMREX_RESTRICT p,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] RefPart const & refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // normalize quad units to MAD-X convention if needed
            amrex::ParticleReal dpx = m_xkick;
            amrex::ParticleReal dpy = m_ykick;
            if (m_unit == UnitSystem::Tm) {
                  dpx /= refpart.rigidity_Tm();
                  dpy /= refpart.rigidity_Tm();
            }

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum
            p.pos(RealAoS::x) = x;
            pxout = px + dpx;

            p.pos(RealAoS::y) = y;
            pyout = py + dpy;

            p.pos(RealAoS::t) = t;
            ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

    private:
        amrex::ParticleReal m_xkick; //! horizontal kick strength
        amrex::ParticleReal m_ykick; //! vertical kick strength
        UnitSystem m_unit; //! Kicks are for 0 dimensionless, or for 1 in T-m."
    };

} // namespace impactx

#endif // IMPACTX_KICKER_H
