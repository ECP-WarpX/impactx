/* Copyright 2021 Axel Huebl, Chad Mitchell
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIAGNOSTIC_OUTPUT_H
#define IMPACTX_DIAGNOSTIC_OUTPUT_H

#include "particles/ImpactXParticleContainer.H"


namespace impactx
{
namespace diagnostics
{
    /** Type of beam diagnostic output
     */
    enum class OutputType
    {
        PrintParticles
    };

    /** Output diagnostics associated with the beam.
     *
     * @param pc container of the particles use for diagnostics
     * @param otype the type of output to produce
     */
    void DiagnosticOutput(ImpactXParticleContainer &pc,
                                  OutputType const & otype);

} // namespace diagnostics
} // namespace impactx

#endif // IMPACTX_DIAGNOSTIC_OUTPUT_H
