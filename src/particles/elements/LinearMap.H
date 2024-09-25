/* Copyright 2022-2024 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENT_LINEAR_MAP_H
#define IMPACTX_ELEMENT_LINEAR_MAP_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/lineartransport.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct LinearMap
    : public elements::BeamOptic<LinearMap>,
      public elements::Thin,
      public elements::Alignment,
      public elements::LinearTransport,
      public elements::NoFinalize
    {
        static constexpr auto type = "LinearMap";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A thin element that applies a user-provided linear transport map R
         *  to the 6-vector of phase space coordinates (x,px,y,py,t,pt).
         *  Thus  x_final = R(1,1)*x + R(1,2)*px + R(1,3)*y + ...,
         *       px_final = R(2,1)*x + R(2,2)*px + R(2,3)*y + ..., etc.
         *
         * @param R user-provided transport map
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        LinearMap (
            LinearTransport::Map6x6 const & R,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0
        )
        : Alignment(dx, dy, rotation_degree)
        {
            for (int i=1; i<=6; ++i) {
                for (int j = 1; j <= 6; ++j) {
                    m_transport_map(i, j) = R(i, j);
                }
            }
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a LinearMap functor, so that a variable of this type can be used like a
         *  LinearMap function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unused)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unused)
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
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
            [[maybe_unused]] RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // input and output phase space vectors
            amrex::Array1D<amrex::ParticleReal,1,6> vectorin;
            amrex::Array1D<amrex::ParticleReal,1,6> vectorout;
            vectorin(1) = x;
            vectorin(2) = px;
            vectorin(3) = y;
            vectorin(4) = py;
            vectorin(5) = t;
            vectorin(6) = pt;

            for (int i=1; i<=6; ++i) {

                vectorout(i) = 0.0;
                for (int j = 1; j <= 6; ++j) {
                    vectorout(i) += m_transport_map(i, j) * vectorin(j);
                }

            }

            // assign updated values
            x = vectorout(1);
            px = vectorout(2);
            y = vectorout(3);
            py = vectorout(4);
            t = vectorout(5);
            pt = vectorout(6);

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        LinearTransport::Map6x6 m_transport_map;  // 6x6 transport map

    };

} // namespace impactx

#endif // IMPACTX_ELEMENT_LINEAR_MAP_H
