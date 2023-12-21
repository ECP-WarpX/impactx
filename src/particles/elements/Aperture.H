/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_APERTURE_H
#define IMPACTX_APERTURE_H


#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"

#include <ablastr/particles/IndexHandling.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace impactx
{
    struct Aperture
    : public elements::BeamOptic<Aperture>,
      public elements::Thin,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto name = "Aperture";
        using PType = ImpactXParticleContainer::ParticleType;

        enum Shape
        {
            rectangular,
            elliptical
        };

        /** A thin collimator element that applies a transverse aperture boundary.
         *  Particles outside the boundary are considered lost.
         *
         * @param shape aperture shape
         * @param xmax maximum value of horizontal coordinate (m)
         * @param ymax maximum value of vertical coordinate (m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error around the s axis (in the x-y plane) [degrees]
         */
        Aperture (
            amrex::ParticleReal xmax,
            amrex::ParticleReal ymax,
            Shape shape,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree),
          m_shape(shape), m_xmax(xmax), m_ymax(ymax)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is an aperture functor, so that a variable of this type can be used like an
         *  aperture function.
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
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT px,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT py,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] RefPart const & refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            auto const id = p.id();

            // scale horizontal and vertical coordinates
            amrex::ParticleReal const u = x / m_xmax;
            amrex::ParticleReal const v = y / m_ymax;

            // compare against the aperture boundary
            switch (m_shape)
            {
                case Shape::rectangular :  // default
                  if (pow(u,2)>1 || pow(v,2) > 1_prt) {
                     p.id() = -id;
                  }
                  break;

               case Shape::elliptical :
                  if (pow(u,2)+pow(v,2) > 1_prt) {
                     p.id() = -id;
                  }
                  break;
            }

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

    private:
        Shape m_shape; //! aperture type (rectangular, elliptical)
        amrex::ParticleReal m_xmax; //! maximum horizontal coordinate
        amrex::ParticleReal m_ymax; //! maximum vertical coordinate

    };

} // namespace impactx

#endif // IMPACTX_APERTURE_H
