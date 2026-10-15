/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Weiqun Zhang
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EMBEDDED_BOUNDARY_H
#define IMPACTX_EMBEDDED_BOUNDARY_H

#include "AmrCoreData_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EBFabFactory.H>
#else
#   include <AMReX_FabFactory.H>
#   include <AMReX_FArrayBox.H>
#endif

#include <array>
#include <optional>
#include <string>


namespace impactx::initialization
{
#ifdef AMREX_USE_EB
    using EBFactory = amrex::EBFArrayBoxFactory;
#else
    using EBFactory = amrex::FArrayBoxFactory;
#endif

    /** A conducting transverse wall (beam pipe), infinite along z */
    struct TransverseWall
    {
        /** implicit function f(x,y): negative inside the pipe, positive in the wall */
        std::string implicit_function;
        /** transverse bounding box (x, y) of the pipe aperture */
        std::array<amrex::Real, 2> bbox_lo;
        std::array<amrex::Real, 2> bbox_hi;
    };

    /** Read the transverse wall from the eb.* input parameters
     *
     * @return the wall, or nullopt if eb.shape is none (the default)
     */
    std::optional<TransverseWall>
    read_transverse_wall ();

    /** Build the EB data of the transverse wall on all levels
     *
     * Needs to be called again whenever the geometry changes.
     * This is a no-op if eb.shape is none.
     */
    void
    build_eb (AmrCoreData & amr_data);

    /** The EB factories per level, or empty if there is no EB */
    amrex::Vector<EBFactory const *>
    eb_factories (AmrCoreData const & amr_data);

    /** Release the EB data on all levels */
    void
    clear_eb (AmrCoreData & amr_data);

} // namespace impactx::initialization

#endif // IMPACTX_EMBEDDED_BOUNDARY_H
