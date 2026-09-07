/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
     *
     * If MPI is enabled and not yet initialized, MPI is initialized here (and
     * finalized at process exit), so that MPI outlives all
     * amrex::Initialize/amrex::Finalize cycles of the process.
     */
    void default_init_AMReX (int argc, char* argv[]);

    /** Initialize AMReX
     *
     * Initializes AMReX if not already done.
     *
     * If MPI is enabled and not yet initialized, MPI is initialized here (and
     * finalized at process exit), so that MPI outlives all
     * amrex::Initialize/amrex::Finalize cycles of the process.
     */
    void default_init_AMReX ();

    /** Finalize AMReX, if ImpactX initialized it
     *
     * This is a no-op if AMReX was initialized by another library or by the
     * user (e.g., via pyAMReX), because then the entity that initialized AMReX
     * also finalizes it.
     *
     * @return true if amrex::Finalize() was called
     */
    bool default_finalize_AMReX ();

} // namespace impactx::initialization

#endif // IMPACT_DEFAULT_INIT_AMREX_H
