/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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

#include "Aperture.H"
#include "Buncher.H"
#include "CFbend.H"
#include "ChrDrift.H"
#include "ChrPlasmaLens.H"
#include "ChrQuad.H"
#include "ChrSoftQuad.H"
#include "ChrUniformAcc.H"
#include "ConstF.H"
#include "DipEdge.H"
#include "Drift.H"
#include "Empty.H"
#include "ExactCFbend.H"
#include "ExactDrift.H"
#include "ExactMultipole.H"
#include "ExactQuad.H"
#include "ExactSbend.H"
#include "Kicker.H"
#include "LinearMap.H"
#include "Marker.H"
#include "Multipole.H"
#include "NonlinearLens.H"
#include "PlaneXYRot.H"
#include "PolygonAperture.H"
#include "Programmable.H"
#include "PRot.H"
#include "Quad.H"
#include "QuadEdge.H"
#include "RFCavity.H"
#include "Sbend.H"
#include "ShortRF.H"
#include "SoftSol.H"
#include "SoftQuad.H"
#include "Sol.H"
#include "SpinMap.H"
#include "Source.H"
#include "TaperedPL.H"
#include "ThinDipole.H"
#include "diagnostics/BeamMonitor.H"

#include <variant>

namespace impactx::elements
{
    using KnownElements = std::variant<
        Empty, /* must be first, so KnownElements creates a default constructor */
        Aperture,
        Buncher,
        CFbend,
        ChrAcc,
        ChrDrift,
        ChrPlasmaLens,
        ChrQuad,
        ChrSoftQuadrupole,
        ConstF,
        diagnostics::BeamMonitor,
        DipEdge,
        Drift,
        ExactCFbend,
        ExactDrift,
        ExactMultipole,
        ExactQuad,
        ExactSbend,
        Kicker,
        LinearMap,
        Marker,
        Multipole,
        NonlinearLens,
        PlaneXYRot,
        PolygonAperture,
        Programmable,
        PRot,
        Quad,
        QuadEdge,
        RFCavity,
        Sbend,
        ShortRF,
        SoftSolenoid,
        SoftQuadrupole,
        Sol,
        Source,
        SpinMap,
        TaperedPL,
        ThinDipole
    >;

} // namespace impactx

#endif // IMPACTX_ELEMENTS_ALL_H
