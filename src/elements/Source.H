/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOURCE_H
#define IMPACTX_SOURCE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/lineartransport.H"
#include "mixin/nofinalize.H"
#include "mixin/thin.H"

#include <optional>
#include <stdexcept>
#include <string>

namespace impactx::elements
{
    struct SourceGPU
    : public mixin::LinearTransport<SourceGPU>,
      public mixin::Thin,
      public mixin::NoFinalize
    {
        static constexpr auto type = "Source";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr bool DEFAULT_active_once = true;
        static constexpr bool DEFAULT_load_ref_particle = true;
        static constexpr std::optional<int> DEFAULT_load_step = std::nullopt;
        static constexpr std::optional<int> DEFAULT_load_step_index = std::nullopt;

        /** A particle source
         *
         * At most one of load_step and load_step_index selects the step to read.
         * If neither is set, the last step in the file is read.
         *
         * @param distribution Must read "openPMD" for this constructor
         * @param openpmd_path path to openPMD series as accepted by openPMD::Series
         * @param active_once inject particles only for the first lattice period
         * @param load_ref_particle restore the reference particle from the species metadata
         *        of the openPMD file (particle tracking only)
         * @param load_step the ImpactX step at which the beam monitor wrote the beam, which is
         *        stored as the openPMD iteration in the file; it must be a step in the file
         * @param load_step_index the position of the step in the file: 0 is the first step and
         *        -1 the last step, counting back from it as in Python (-2 is the second to last)
         * @throws std::invalid_argument if distribution is not "openPMD", if both load_step and
         *         load_step_index are set or if load_step is negative
         */
        SourceGPU (
            std::string distribution,
            std::string openpmd_path,
            bool active_once = DEFAULT_active_once,
            bool load_ref_particle = DEFAULT_load_ref_particle,
            std::optional<int> load_step = DEFAULT_load_step,
            std::optional<int> load_step_index = DEFAULT_load_step_index
        )
        : m_distribution(distribution),
          m_series_name(std::move(openpmd_path)),
          m_active_once(active_once),
          m_load_ref_particle(load_ref_particle),
          m_load_step(load_step),
          m_load_step_index(load_step_index)
        {
            if (distribution != "openPMD") {
                throw std::invalid_argument("Only 'openPMD' distribution is supported if openpmd_path is provided!");
            }

            check_step_selection();
        }

        /** Check that the step to read from the file is selected unambiguously.
         *
         * This is checked on construction and again when particles are loaded,
         * because both options can also be set after the element was created.
         *
         * @throws std::invalid_argument if both load_step and load_step_index are
         *         set or if load_step is negative
         */
        void check_step_selection () const;

        /** Reverse the element (not supported for Source) */
        void reverse ()
        {
            throw std::runtime_error("Source: reverse() is not supported.");
        }

        /** Initialize/Load particles.
         *
         * Particles are relative to the reference particle.
         * The step that is read from the file is selected with load_step or load_step_index.
         * If load_ref_particle is set, the reference particle is restored from
         * the species metadata of the file before the particles are added.
         *
         * @param[in,out] pc particle container to push
         * @param[in] step current global step of the tracking loop (unused)
         * @param[in] period for periodic lattices, this is the current period (turn or cycle)
         * @throws std::invalid_argument if both load_step and load_step_index are set or if
         *         load_step is negative
         * @throws std::runtime_error if the series holds no step or if the selected step or
         *         index is not in the file
         */
        void operator() (
            ImpactXParticleContainer & pc,
            [[maybe_unused]] int step,
            int period
        );

        /** This function returns the linear transport map.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            // nothing to do
            return Map6x6::Identity();
        }

        /** This does nothing to the reference particle. */
        using Thin::operator();

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        std::string m_distribution; //! Distribution type of particles in the source
        std::string m_series_name; //! openPMD filename
        bool m_active_once = true; //! Inject particles only for the first lattice period
        bool m_load_ref_particle = true; //! Restore the reference particle from the file's species metadata
        std::optional<int> m_load_step = DEFAULT_load_step; //! ImpactX step (openPMD iteration) to load
        std::optional<int> m_load_step_index = DEFAULT_load_step_index; //! position of the step in the file, -1 is the last step
    };

    /** Source: GPU state plus common host-only metadata */
    using Source = Element<SourceGPU>;

} // namespace impactx

#endif  // IMPACTX_SOURCE_H
