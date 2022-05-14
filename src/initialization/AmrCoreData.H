/* Copyright 2021-2022 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_AMR_CORE_DATA_H
#define IMPACTX_AMR_CORE_DATA_H

#include <AMReX_AmrCore.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>


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

        AmrCoreData (AmrCoreData&& rhs) = default;
        AmrCoreData& operator= (AmrCoreData&& rhs) = default;

        AmrCoreData (const AmrCoreData& rhs) = delete;
        AmrCoreData& operator= (const AmrCoreData& rhs) = delete;

        ~AmrCoreData() = default;

    private:
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
