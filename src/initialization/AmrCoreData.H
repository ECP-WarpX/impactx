/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_AMR_CORE_DATA_H
#define IMPACTX_AMR_CORE_DATA_H

#include "AmrCoreData_fwd.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/CovarianceMatrix.H"

#include <AMReX_AmrCore.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EB2.H>
#   include <AMReX_EBFabFactory.H>
#endif
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>

#include <memory>
#include <optional>
#include <string>
#include <unordered_map>


namespace impactx::initialization
{
    /** Data in AmrCore
     *
     * This wrapper is needed because AmrCore cannot be constructed, due to
     * pure virtuals in it.
     *
     * The primary purpose of this wrapper is to transport simulation_geometry
     * the geometry (topology) of the simulation and amr_info; contains
     * information on mesh-refinement and box/grid blocks.
     */
    class AmrCoreData final
        : public amrex::AmrCore
    {
      public:
        AmrCoreData (
            amrex::Geometry const& level_0_geom,
            amrex::AmrInfo const& amr_info
        );
        AmrCoreData (
            amrex::RealBox const & rb,
            int max_level_in,
            amrex::Vector<int> const & n_cell_in,
            int coord,
            amrex::Vector<amrex::IntVect> const & ref_ratios,
            amrex::Array<int,AMREX_SPACEDIM> const & is_per
        );

        AmrCoreData (AmrCoreData&& rhs) = default;
        AmrCoreData& operator= (AmrCoreData&& rhs) = default;

        AmrCoreData (const AmrCoreData& rhs) = delete;
        AmrCoreData& operator= (const AmrCoreData& rhs) = delete;

        ~AmrCoreData() override;

        /** Data needed for the symplectic particle tracking mode */
        struct TrackModeParticles
        {
            /** these are the physical/beam particles of the simulation */
            std::unique_ptr<impactx::ImpactXParticleContainer> m_particle_container;

            /** former beam particles that got lost in apertures, the wall, etc. */
            std::unique_ptr<impactx::ImpactXParticleContainer> m_particles_lost;

            /** charge density per level */
            std::unordered_map<int, amrex::MultiFab> m_rho;
            /** charge density projected to 2D (x-y) per level, for the 2D
             *  space-charge models; (re)derived from m_rho on access */
            std::unordered_map<int, amrex::MultiFab> m_rho_2d;
            /** scalar potential per level */
            std::unordered_map<int, amrex::MultiFab> m_phi;
            /** space charge field (vector) per level */
            std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > m_space_charge_field;

#ifdef AMREX_USE_EB
            /** EB factory of the transverse wall per level */
            std::unordered_map<int, std::unique_ptr<amrex::EBFArrayBoxFactory> > m_eb_factory;
            /** EB index space per level, owned by the global EB2 stack */
            std::unordered_map<int, amrex::EB2::IndexSpace*> m_eb_index_space;
#endif

        } track_particles;

        /** Data needed for the beam envelop (covariance matrix) tracking mode */
        struct TrackModeEnvelope
        {
            /** The reference particle */
            std::optional<RefPart> m_ref;

            /** The 6x6 covariance matrix for the beam envelope, relative to the reference particle */
            std::optional<Envelope> m_env;

        } track_envelope;

        /** Data needed for the reference orbit tracking mode */
        struct TrackModeReference
        {
            /** The reference particle */
            std::optional<RefPart> m_ref;

        } track_reference;

        void ErrorEst (
            [[maybe_unused]] int lev,
            [[maybe_unused]] amrex::TagBoxArray& tags,
            [[maybe_unused]] amrex::Real time,
            [[maybe_unused]] int ngrow) override;

        void MakeNewLevelFromScratch (
            [[maybe_unused]] int lev,
            [[maybe_unused]] amrex::Real time,
            [[maybe_unused]] const amrex::BoxArray& ba,
            [[maybe_unused]] const amrex::DistributionMapping& dm) override;

        void MakeNewLevelFromCoarse (
            [[maybe_unused]] int lev,
            [[maybe_unused]] amrex::Real time,
            [[maybe_unused]] const amrex::BoxArray& ba,
            [[maybe_unused]] const amrex::DistributionMapping& dm) override;

        void RemakeLevel (
            [[maybe_unused]] int lev,
            [[maybe_unused]] amrex::Real time,
            [[maybe_unused]] const amrex::BoxArray& ba,
            [[maybe_unused]] const amrex::DistributionMapping& dm) override;

        void ClearLevel ([[maybe_unused]] int lev) override;
    };

} // namespace impactx::initialization

#endif // IMPACTX_AMR_CORE_DATA_H
