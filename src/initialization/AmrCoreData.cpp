/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "AmrCoreData.H"

#include <AMReX.H>


namespace impactx::initialization
{
    AmrCoreData::AmrCoreData (
        amrex::Geometry const& level_0_geom,
        amrex::AmrInfo const& amr_info
    )
        : amrex::AmrCore(level_0_geom, amr_info)
    {
    }

    AmrCoreData::AmrCoreData (
        amrex::RealBox const & rb,
        int max_level_in,
        amrex::Vector<int> const & n_cell_in,
        int coord,
        amrex::Vector<amrex::IntVect> const & ref_ratios,
        amrex::Array<int,AMREX_SPACEDIM> const & is_per
    )
        : amrex::AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per)
    {
    }

    void
    AmrCoreData::ErrorEst (
        [[maybe_unused]] int lev,
        [[maybe_unused]] amrex::TagBoxArray& tags,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] int ngrow)
    {
        amrex::Abort("Do not use");
    }

    void
    AmrCoreData::MakeNewLevelFromScratch (
        [[maybe_unused]] int lev,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] const amrex::BoxArray& ba,
        [[maybe_unused]] const amrex::DistributionMapping& dm)
    {
        amrex::Abort("Do not use");
    }

    void
    AmrCoreData::MakeNewLevelFromCoarse (
        [[maybe_unused]] int lev,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] const amrex::BoxArray& ba,
        [[maybe_unused]] const amrex::DistributionMapping& dm)
    {
        amrex::Abort("Do not use");
    }

    void
    AmrCoreData::RemakeLevel (
        [[maybe_unused]] int lev,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] const amrex::BoxArray& ba,
        [[maybe_unused]] const amrex::DistributionMapping& dm)
    {
        amrex::Abort("Do not use");
    }

    void
    AmrCoreData::ClearLevel ([[maybe_unused]] int lev)
    {
        amrex::Abort("Do not use");
    }
} // namespace impactx::initialization
