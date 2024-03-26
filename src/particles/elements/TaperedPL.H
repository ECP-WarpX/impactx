/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TAPEREDPL_H
#define IMPACTX_TAPEREDPL_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct TaperedPL
    : public elements::BeamOptic<TaperedPL>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "TaperedPL";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A short segment of a nonlinear plasma lens with a transverse taper.
         *
         * @param k  integrated plasma lens (linear) focusing strength in m^(-1)
         *           = (length in m) * (azimuthal field gradient in T/m) / (rigidity in T-m)
         *        OR integrated azimuthal magnetic field gradient in T (k > 0)
         * @param unit  Unit specification
         *           unit = 0 integrated linear focusing strength in m^(-1)
         *           unit = 1 integrated linear focusing strength in T
         * @param taper  horizontal taper parameter in m^(-1)
                    = 1 / (horizontal dispersion in m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        TaperedPL (
            amrex::ParticleReal k,
            amrex::ParticleReal taper,
            int unit,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_k(k), m_taper(taper), m_unit(unit)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a tapered PL functor, so that a variable of this type can be used like a
         *  shortrf function.
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

            // normalize focusing strength units to MAD-X convention if needed
            amrex::ParticleReal g = m_k;
            if (m_unit == 1) {
                  g = m_k / refpart.rigidity_Tm();
            }

            // intialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum
            // xout = x;
            pxout = px - g * ( x + m_taper/2.0 * (pow(x,2) + pow(y,2)) );

            // yout = y;
            pyout = py - g * ( y + m_taper * x * y );

            // tout = t;
            ptout = pt;

            // assign updated values
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        amrex::ParticleReal m_k; //! linear focusing strength (field gradient)
        amrex::ParticleReal m_taper; //! horizontal taper parameter
        int m_unit;  //! units for linear focusing strength (field gradient)
    };

} // namespace impactx

#endif // IMPACTX_TAPEREDPL_H
