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

#include <string>


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

    /** Read the parameters of an inputs file
     *
     * Initializes AMReX if not already done. In that case, the parameters of the file are
     * read before AMReX and ImpactX look up their runtime parameters, as for an inputs file
     * passed to the impactx executable. If that fails, AMReX is not initialized.
     *
     * @param filename path to the inputs file
     * @throws std::filesystem::filesystem_error if the file cannot be read
     * @throws std::runtime_error if AMReX finds an error in the file
     */
    void load_inputs_file (std::string const & filename);

} // namespace impactx::initialization

#endif // IMPACT_DEFAULT_INIT_AMREX_H
