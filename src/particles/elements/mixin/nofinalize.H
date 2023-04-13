/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_NOFINALIZE_H
#define IMPACTX_ELEMENTS_MIXIN_NOFINALIZE_H


namespace impactx::elements
{
    /** This is a helper class for lattice elements that need no finalize function.
     *
     * Finalize helpers are used to clean up static/shared data, usually dynamic
     * memory from parameters or state from external libraries.
     */
    struct NoFinalize
    {
        /** Close and deallocate all data and handles.
         */
        void
        finalize ()
        {
            // nothing to do
        }
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_NOFINALIZE_H
