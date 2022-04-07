/* Copyright 2021-2022 Axel Huebl
 *
 * This file is part of ImpactX.
 *
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
