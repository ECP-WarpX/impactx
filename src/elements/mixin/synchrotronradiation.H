/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_SYNCHROTRON_RADIATION_H
#define IMPACTX_ELEMENTS_MIXIN_SYNCHROTRON_RADIATION_H


namespace impactx::elements::mixin
{
    /** This is a helper class for elements that emit synchrotron radiation
     *
     * An element deriving from this class bends the reference orbit, so coherent (CSR) and
     * incoherent (ISR) synchrotron radiation are applied in it when enabled.
     *
     * If an element derives from this class, then it must implement ``rc(refpart)``, the
     * signed radius of curvature of the reference orbit in m.
     */
    struct SynchrotronRadiation
    {
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_SYNCHROTRON_RADIATION_H
