/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_INIT_PARSER_H
#define IMPACT_INIT_PARSER_H


namespace impactx::initialization
{
    /** Overwrite defaults in AMReX Inputs
     *
     * This overwrites defaults in amrex::ParamParse for inputs.
     */
    void
    overwrite_amrex_parser_defaults ();

} // namespace impactx

#endif // IMPACT_INIT_PARSER_H
