/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_ALL_H
#define IMPACTX_ELEMENTS_ALL_H

#include "ConstF.H"
#include "DipEdge.H"
#include "Drift.H"
#include "Multipole.H"
#include "None.H"
#include "NonlinearLens.H"
#include "Programmable.H"
#include "Quad.H"
#include "RFCavity.H"
#include "Sbend.H"
#include "ShortRF.H"
#include "Sol.H"
#include "PRot.H"
#include "SoftSol.H"
#include "SoftQuad.H"
#include "ChrQuad.H"
#include "ChrDrift.H"
#include "diagnostics/openPMD.H"

#include <variant>


namespace impactx
{
    using KnownElements = std::variant<
        None, /* must be first, so KnownElements creates a default constructor */
        ConstF, DipEdge, Drift, Multipole, NonlinearLens,
        diagnostics::BeamMonitor, Programmable, ChrQuad, ChrDrift,
        Quad, RFCavity, Sbend, ShortRF, Sol, PRot, SoftSolenoid, SoftQuadrupole>;

} // namespace impactx

#endif // IMPACTX_ELEMENTS_ALL_H
