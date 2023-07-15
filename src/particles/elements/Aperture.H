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
      public elements::NoFinalize
    {
        static constexpr auto name = "Aperture";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A thin collimator element that applies a transverse aperture boundary.
         *  Particles outside the boundary are considered lost.
         *
         * @param shape aperture shape (0 rectangular, 1 elliptical)
         * @param xmax maximum value of horizontal coordinate (m)
         * @param ymax maximum value of vertical coordinate (m)
         */
        Aperture( int const shape,
                   amrex::ParticleReal const xmax,
                   amrex::ParticleReal const ymax )
        : m_shape(shape), m_xmax(xmax), m_ymax(ymax)
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
                [[maybe_unused]] RefPart const & refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            int const id = p.id();
            //uint64_t const global_id = ablastr::particles::localIDtoGlobal(p.id(), p.cpu());

            // scale horizontal and vertical coordinates
            amrex::ParticleReal const u = x/m_xmax;
            amrex::ParticleReal const v = y/m_ymax;

            // compare against the aperture boundary

            switch(m_shape) {

           // rectangular aperture (default)
               case 0 :
                  if (pow(u,2)>1 || pow(v,2)>1){
                     p.pos(RealAoS::x) = 0.0_prt; //replace with id change of sign
                     p.pos(RealAoS::y) = 0.0_prt; //replace with id change of sign
                     p.id() = -id;
                  }
                  break;

           // elliptical aperture
               case 1 :
                  if (pow(u,2)+pow(v,2)>1){
                     p.pos(RealAoS::x) = 0.0_prt; //replace with id change of sign
                     p.pos(RealAoS::y) = 0.0_prt; //replace with id change of sign
                     p.id() = -id;
                  }
                  break;

            }

        }

        /** This pushes the reference particle. */
        using Thin::operator();

    private:
        int m_shape; //! aperture type (1 rectangular, 2 elliptical)
        amrex::ParticleReal m_xmax; //! maximum horizontal coordinate
        amrex::ParticleReal m_ymax; //! maximum vertical coordinate

    };

} // namespace impactx

#endif // IMPACTX_APERTURE_H
