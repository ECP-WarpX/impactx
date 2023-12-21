/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_ALIGNMENT_H
#define IMPACTX_ELEMENTS_MIXIN_ALIGNMENT_H

#include "particles/ImpactXParticleContainer.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx::elements
{
    /** This is a helper class for lattice elements with horizontal/vertical alignment errors
     */
    struct Alignment
    {
        static constexpr amrex::ParticleReal degree2rad = ablastr::constant::math::pi / 180.0;

        /** A finite-length element
         *
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error around the s axis (in the x-y plane) [degrees]
         */
        Alignment (
            amrex::ParticleReal dx,
            amrex::ParticleReal dy,
            amrex::ParticleReal rotation_degree
        )
        : m_dx(dx), m_dy(dy), m_rotation(rotation_degree * degree2rad)
        {
        }

        Alignment () = default;
        Alignment (Alignment const &) = default;
        Alignment& operator= (Alignment const &) = default;
        Alignment (Alignment&&) = default;
        Alignment& operator= (Alignment&& rhs) = default;

        ~Alignment () = default;

        /** Shift the particle into the alignment error frame
         *
         * @param[inout] x horizontal position relative to reference particle
         * @param[inout] y vertical position relative to reference particle
         * @param[inout] px horizontal momentum relative to reference particle
         * @param[inout] py vertical momentum relative to reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void shift_in (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py
        ) const
        {
            auto const [sin_rotation, cos_rotation] = amrex::Math::sincos(m_rotation);

            // position
            amrex::ParticleReal const xc = x - m_dx;
            amrex::ParticleReal const yc = y - m_dy;
            x =  xc * cos_rotation + yc * sin_rotation;
            y = -xc * sin_rotation + yc * cos_rotation;

            // momentum
            amrex::ParticleReal const pxc = px;
            amrex::ParticleReal const pyc = py;
            px =  pxc * cos_rotation + pyc * sin_rotation;
            py = -pxc * sin_rotation + pyc * cos_rotation;
        }

        /** Shift the particle out of the alignment error frame
         *
         * @param[inout] x horizontal position relative to reference particle
         * @param[inout] y vertical position relative to reference particle
         * @param[inout] px horizontal momentum relative to reference particle
         * @param[inout] py vertical momentum relative to reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void shift_out (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py
        ) const
        {
            auto const [sin_rotation, cos_rotation] = amrex::Math::sincos(m_rotation);

            // position
            amrex::ParticleReal const xc = x;
            amrex::ParticleReal const yc = y;
            x = xc * cos_rotation - yc * sin_rotation;
            y = xc * sin_rotation + yc * cos_rotation;
            x += m_dx;
            y += m_dy;

            // momentum
            amrex::ParticleReal const pxc = px;
            amrex::ParticleReal const pyc = py;
            px = pxc * cos_rotation - pyc * sin_rotation;
            py = pxc * sin_rotation + pyc * cos_rotation;
        }

        /** Horizontal translation error
         *
         * @return horizontal translation error in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal dx () const
        {
            return m_dx;
        }

        /** Vertical translation error
         *
         * @return vertical translation error in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal dy () const
        {
            return m_dy;
        }

        /** rotation error around the s axis (in the x-y plane)
         *
         * @return rotation error around the s axis (in the x-y plane) [degrees]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal rotation () const
        {
            return m_rotation / degree2rad;
        }

        amrex::ParticleReal m_dx = 0; //! horizontal translation error [m]
        amrex::ParticleReal m_dy = 0; //! vertical translation error [m]
        amrex::ParticleReal m_rotation = 0; //! rotation error around the s axis (in the x-y plane) [rad]
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_ALIGNMENT_H
