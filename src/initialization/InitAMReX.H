/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_DEFAULT_INIT_AMREX_H
#define IMPACT_DEFAULT_INIT_AMREX_H

#include "AmrCoreData.H"


namespace impactx::initialization
{
    /** Initialize AMReX
     *
     * Initializes AMReX if not already done.
     */
    void default_init_AMReX (int argc, char* argv[]);

    /** Initialize AMReX
     *
     * Initializes AMReX if not already done.
     */
    void default_init_AMReX ();

} // namespace impactx::initialization

#endif // IMPACT_DEFAULT_INIT_AMREX_H
