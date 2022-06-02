/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_INIT_ONE_BOX_PER_RANK_H
#define IMPACT_INIT_ONE_BOX_PER_RANK_H

#include "AmrCoreData.H"


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
    AmrCoreData
    one_box_per_rank ();

} // namespace impactx::initialization

#endif // IMPACT_INIT_ONE_BOX_PER_RANK_H
