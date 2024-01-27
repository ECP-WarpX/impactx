/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_AMR_CORE_DATA_FWD_H
#define IMPACTX_AMR_CORE_DATA_FWD_H


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
    class AmrCoreData;

} // namespace impactx::initialization

#endif // IMPACTX_AMR_CORE_DATA_FWD_H
