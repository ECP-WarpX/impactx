/* Copyright 2021-2022 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_INIT_ONE_BOX_PER_RANK_H
#define IMPACT_INIT_ONE_BOX_PER_RANK_H

#include <AMReX_Geometry.H>
#include <AMReX_AmrCore.H>

#include <tuple>


namespace impactx::initialization
{
    /** This builds one AMReX box per MPI-rank
     *
     * This is a simple decomposition for particles that we default to if
     * we do not need space charge or load balancing.
     *
     * @return simulation_geometry the geometry (topology) of the simulation;
     *         amr_info contains information on mesh-refinement and box/grid blocks
     */
    //std::tuple< amrex::Geometry, amrex::AmrInfo >
    amrex::AmrMesh
    one_box_per_rank ();

} // namespace impactx::initialization

#endif // IMPACT_INIT_ONE_BOX_PER_RANK_H
