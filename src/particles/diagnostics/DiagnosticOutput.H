/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIAGNOSTIC_OUTPUT_H
#define IMPACTX_DIAGNOSTIC_OUTPUT_H

#include "particles/ImpactXParticleContainer.H"

#include <string>


namespace impactx::diagnostics
{
    /** Type of beam diagnostic output
     */
    enum class OutputType
    {
        PrintParticles, ///< ASCII diagnostics, for small tests only
        PrintNonlinearLensInvariants, ///< ASCII diagnostics for the IOTA nonlinear lens, for small tests only
        PrintRefParticle, ///< ASCII diagnostics, for small tests only
        PrintReducedBeamCharacteristics ///< ASCII diagnostics, for beam momenta and Twiss parameters
    };

    /** ASCII output diagnostics associated with the beam.
     *
     * This temporary implementation uses ASCII output.
     * It is intended only for small tests where IO performance is not
     * a concern. The implementation here serializes IO.
     *
     * @param pc container of the particles use for diagnostics
     * @param otype the type of output to produce
     * @param file_name the file name to write to
     * @param step the global step
     * @param append open a new file with a fresh header (false) or append data to an existing file (true)
     */
    void DiagnosticOutput (ImpactXParticleContainer const & pc,
                           OutputType const otype,
                           std::string file_name,
                           int const step = 0,
                           bool const append = false);

} // namespace impactx::diagnostics

#endif // IMPACTX_DIAGNOSTIC_OUTPUT_H
