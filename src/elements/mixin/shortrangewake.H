/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_SHORT_RANGE_WAKE_H
#define IMPACTX_ELEMENTS_MIXIN_SHORT_RANGE_WAKE_H

#include <AMReX_REAL.H>
#include <AMReX_TrackedVector.H>

#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>


namespace impactx::elements::mixin
{
    /** Non-owning view on the user-specified single-particle wake function of an element.
     *
     * Trivially copyable, so it can be captured by value into GPU kernels. The
     * pointers address memory owned by the element and are only valid for the
     * duration of the current call: do not cache them, re-query per use instead.
     */
    struct ShortRangeWakeView
    {
        int n_data = 0;  //!< number of tabulated values per wake function
        amrex::ParticleReal dz = 0;  //!< spacing of the longitudinal positions in m
        amrex::ParticleReal const * z = nullptr;  //!< longitudinal positions in m
        amrex::ParticleReal const * wake_x = nullptr;  //!< horizontal single-particle wake function
        amrex::ParticleReal const * wake_y = nullptr;  //!< vertical single-particle wake function
        amrex::ParticleReal const * wake_z = nullptr;  //!< longitudinal single-particle wake function
    };

    /** This is a helper class for elements with a short-range structure wake
     *
     * It holds the tabulated single-particle wake function the user specified for the
     * element. Host data persists across simulations; device data is uploaded on-demand
     * "lazily" on first use and released when AMReX shuts down.
     *
     * The wake is optional: an element without tabulated data does not contribute a
     * short-range wake. An element deriving from this class calls @see release_wake_gpu
     * from its ``finalize()``.
     */
    struct ShortRangeWake
    {
        ShortRangeWake () = default;

        /** An element with an optional user-specified wake function
         *
         * Absent data is stored as empty. With z_data, wake_z_data has one value per position.
         *
         * @param z_data longitudinal positions of the single-particle wake function in m,
         *               at least two
         * @param wake_x_data horizontal single-particle wake function data
         * @param wake_y_data vertical single-particle wake function data
         * @param wake_z_data longitudinal single-particle wake function data
         * @throws std::invalid_argument if z_data has fewer than two values, its first two
         *         values are equal, or wake_z_data has a different number of values
         */
        ShortRangeWake (
            std::optional<std::vector<amrex::ParticleReal>> z_data,
            std::optional<std::vector<amrex::ParticleReal>> wake_x_data,
            std::optional<std::vector<amrex::ParticleReal>> wake_y_data,
            std::optional<std::vector<amrex::ParticleReal>> wake_z_data
        )
          : m_z_data(std::move(z_data).value_or(std::vector<amrex::ParticleReal>{})),
            m_wake_x_data(std::move(wake_x_data).value_or(std::vector<amrex::ParticleReal>{})),
            m_wake_y_data(std::move(wake_y_data).value_or(std::vector<amrex::ParticleReal>{})),
            m_wake_z_data(std::move(wake_z_data).value_or(std::vector<amrex::ParticleReal>{}))
        {
            auto const & z = m_z_data.host_const();
            if (z.empty()) { return; }

            if (z.size() < 2) {
                throw std::invalid_argument("ShortRangeWake: z_data needs at least two values!");
            }
            if (m_wake_z_data.host_const().size() != z.size()) {
                throw std::invalid_argument(
                    "ShortRangeWake: wake_z_data needs as many values as z_data!");
            }

            // the spacing of the tabulated positions, from the host data: the push only gets
            // device pointers
            m_dz = std::abs(z[1] - z[0]);
            if (m_dz == 0) {
                throw std::invalid_argument(
                    "ShortRangeWake: the first two values of z_data must differ!");
            }
        }

        /** Does this element carry a tabulated wake function? */
        [[nodiscard]] bool has_short_range_wake () const
        {
            return !m_z_data.host_const().empty();
        }

        /** A view on the wake function in the execution space of the push
         *
         * Uploads the data to the device on first use (a no-op on CPU builds).
         */
        [[nodiscard]] ShortRangeWakeView short_range_wake () const
        {
            ShortRangeWakeView view;
            auto const & d_z = m_z_data.device_const();
            view.n_data = static_cast<int>(d_z.size());
            view.dz = m_dz;
            view.z = d_z.data();
            view.wake_x = m_wake_x_data.device_const().data();
            view.wake_y = m_wake_y_data.device_const().data();
            view.wake_z = m_wake_z_data.device_const().data();
            return view;
        }

        /** Longitudinal positions of the wake function in m */
        [[nodiscard]] std::vector<amrex::ParticleReal> z_data () const { return m_z_data.host_const(); }
        /** Horizontal single-particle wake function */
        [[nodiscard]] std::vector<amrex::ParticleReal> wake_x_data () const { return m_wake_x_data.host_const(); }
        /** Vertical single-particle wake function */
        [[nodiscard]] std::vector<amrex::ParticleReal> wake_y_data () const { return m_wake_y_data.host_const(); }
        /** Longitudinal single-particle wake function */
        [[nodiscard]] std::vector<amrex::ParticleReal> wake_z_data () const { return m_wake_z_data.host_const(); }

        /** Release device storage; host data survives for reuse */
        void release_wake_gpu ()
        {
            m_z_data.release_gpu();
            m_wake_x_data.release_gpu();
            m_wake_y_data.release_gpu();
            m_wake_z_data.release_gpu();
        }

    private:
        amrex::ParticleReal m_dz = 0;  //!< spacing of the longitudinal positions in m
        amrex::Gpu::TrackedVector<amrex::ParticleReal> m_z_data;  //!< longitudinal positions in m
        amrex::Gpu::TrackedVector<amrex::ParticleReal> m_wake_x_data;  //!< horizontal wake function
        amrex::Gpu::TrackedVector<amrex::ParticleReal> m_wake_y_data;  //!< vertical wake function
        amrex::Gpu::TrackedVector<amrex::ParticleReal> m_wake_z_data;  //!< longitudinal wake function
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_SHORT_RANGE_WAKE_H
