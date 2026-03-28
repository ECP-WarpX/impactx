/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_DEVICEBACKEDARRAYS_H
#define IMPACTX_ELEMENTS_MIXIN_DEVICEBACKEDARRAYS_H

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_TrackedVector.H>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace impactx::elements::mixin
{
    /** A non-owning, device-copyable view of two parallel element arrays.
     *
     * This is what a particle-push kernel sees. It is a plain aggregate of two pointers
     * and a length, so it is cheap to copy into a launch closure and carries no ownership
     * whatsoever. Its validity is bounded by the launch it was prepared for.
     */
    struct DualArrayView
    {
        amrex::ParticleReal const * first = nullptr;
        amrex::ParticleReal const * second = nullptr;
        int size = 0;
    };

    /** Two equal-length arrays of element data, owned directly by the element.
     *
     * Every ImpactX element that carries dynamic data carries exactly this shape: a pair
     * of equal-length ``amrex::ParticleReal`` arrays, whether those are Fourier cosine and
     * sine coefficients, normal and skew multipole coefficients, or polygon x and y
     * vertices. One owner therefore serves all of them, and no per-element lifecycle code
     * is needed.
     *
     * Ownership is ordinary RAII:
     *
     * - a value copy gets independent host arrays, because ``amrex::Gpu::TrackedVector``
     *   deep-copies rather than sharing;
     * - device storage is materialized lazily, on first use inside an AMReX session, and
     *   released on AMReX shutdown while the host arrays survive for reuse;
     * - destruction, explicit release and mutation each wait for earlier device consumers
     *   before touching the arrays.
     */
    class DeviceBackedArrays
    {
    public:
        DeviceBackedArrays () = default;

        /** Own a pair of equal-length arrays
         *
         * @param first the first array
         * @param second the second array, which must have the same length
         * @param what element and array names, used in the error message
         */
        DeviceBackedArrays (
            std::vector<amrex::ParticleReal> first,
            std::vector<amrex::ParticleReal> second,
            std::string const & what
        )
        {
            validate(first, second, what);
            m_first = amrex::Gpu::TrackedVector<amrex::ParticleReal>(std::move(first));
            m_second = amrex::Gpu::TrackedVector<amrex::ParticleReal>(std::move(second));
        }

        DeviceBackedArrays (DeviceBackedArrays const &) = default;
        DeviceBackedArrays (DeviceBackedArrays &&) noexcept = default;

        DeviceBackedArrays & operator= (DeviceBackedArrays const & other)
        {
            if (this != &other)
            {
                wait_for_device_consumers();
                m_first = other.m_first;
                m_second = other.m_second;
                m_device_materialized = false;
            }
            return *this;
        }

        DeviceBackedArrays & operator= (DeviceBackedArrays && other) noexcept
        {
            if (this != &other)
            {
                wait_for_device_consumers();
                m_first = std::move(other.m_first);
                m_second = std::move(other.m_second);
                m_device_materialized = std::exchange(other.m_device_materialized, false);
            }
            return *this;
        }

        ~DeviceBackedArrays ()
        {
            // A kernel that was handed a view of these arrays may still be in flight.
            wait_for_device_consumers();
        }

        /** Common length of both arrays */
        [[nodiscard]] int size () const noexcept
        {
            return static_cast<int>(m_first.host_const().size());
        }

        /** A view of the host arrays, for use in host code
         *
         * Obtain this at the point of use and do not cache it: the arrays may be replaced
         * by a later mutation. Host consumers are the reference-particle push, transfer
         * maps, getters and serialization, none of which run inside a kernel.
         */
        [[nodiscard]] DualArrayView host_view () const
        {
            auto const & first = m_first.host_const();
            auto const & second = m_second.host_const();
            return {first.data(), second.data(), static_cast<int>(first.size())};
        }

        /** A view of the execution-space arrays, for a bounded launch lifetime
         *
         * On a GPU build this materializes the device mirror on first use. On a CPU build
         * it resolves to the host arrays.
         */
        [[nodiscard]] DualArrayView device_view () const
        {
            auto const & first = m_first.device_const();
            auto const & second = m_second.device_const();
            m_device_materialized = true;
            return {first.data(), second.data(), static_cast<int>(first.size())};
        }

        /** A copy of the first array's host values, e.g. for serialization */
        [[nodiscard]] std::vector<amrex::ParticleReal> first_host () const
        {
            return m_first.host_const();
        }

        /** A copy of the second array's host values, e.g. for serialization */
        [[nodiscard]] std::vector<amrex::ParticleReal> second_host () const
        {
            return m_second.host_const();
        }

        /** Replace both arrays as one transaction
         *
         * Values are validated before anything is modified, so a rejected update leaves
         * the element untouched rather than half-applied.
         *
         * @param first the new first array
         * @param second the new second array, which must have the same length
         * @param what element and array names, used in the error message
         */
        void commit (
            std::vector<amrex::ParticleReal> first,
            std::vector<amrex::ParticleReal> second,
            std::string const & what
        )
        {
            validate(first, second, what);
            wait_for_device_consumers();

            m_first = amrex::Gpu::TrackedVector<amrex::ParticleReal>(std::move(first));
            m_second = amrex::Gpu::TrackedVector<amrex::ParticleReal>(std::move(second));
            m_device_materialized = false;
        }

        /** Release device storage, keeping the host arrays for reuse
         *
         * Idempotent. This is what an element's ``finalize()`` delegates to, so that an
         * element survives an AMReX finalize/initialize cycle with its configuration
         * intact and simply re-uploads on next use.
         */
        void release_gpu ()
        {
            wait_for_device_consumers();
            m_first.release_gpu();
            m_second.release_gpu();
            m_device_materialized = false;
        }

    private:
        static void validate (
            std::vector<amrex::ParticleReal> const & first,
            std::vector<amrex::ParticleReal> const & second,
            std::string const & what
        )
        {
            if (first.size() != second.size())
            {
                throw std::invalid_argument(what + " arrays must have the same length!");
            }
        }

        /** Wait for kernels that may still hold a device view of these arrays
         *
         * All streams are synchronized, because different tiles may have been launched on
         * different streams. This is a control-plane cost on destruction, release and
         * mutation only; it never runs during ordinary read-only traversal, and
         * correctness does not depend on anyone else having synchronized first.
         */
        void wait_for_device_consumers () const noexcept
        {
#ifdef AMREX_USE_GPU
            if (m_device_materialized && amrex::Initialized())
            {
                amrex::Gpu::streamSynchronizeAll();
            }
#else
            // no device mirror exists, so nothing can be in flight
#endif
        }

        amrex::Gpu::TrackedVector<amrex::ParticleReal> m_first;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> m_second;

        //! mutable: device_view() is logically const but materializes the mirror
        mutable bool m_device_materialized = false;
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_DEVICEBACKEDARRAYS_H
