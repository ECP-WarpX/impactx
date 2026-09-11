/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TRACKING_STATE_H
#define IMPACTX_TRACKING_STATE_H

#include "elements/helper/Handle.H"
#include "tracking/TrackingDirection.H"

#include <optional>


namespace impactx
{
    /** This class wraps a few states during tracking loops.
     *
     * These states are useful to know in hooks (Python callback functions).
     */
    struct TrackingState
    {
        /** For tracking hooks/callbacks, a global step of the simulation
         *
         * A state of internal simulation steps, increments also for space charge slice steps in elements.
         * We start in "step 0" (initial state).
         */
        int m_step = 0;

        /** For tracking hooks/callbacks, the period in the lattice (e.g., turn or channel period) */
        int m_period = 0;

        /** For tracking hooks/callbacks, the current lattice element
         *
         * A share of the element, not a pointer into the lattice: a hook may hold on to it,
         * and the lattice may be edited afterwards, without leaving a dangling reference.
         */
        std::optional<elements::ElementHandle> m_element;

        /** Tracking direction in the lattice */
        TrackingDirection m_direction = TrackingDirection::Forward;

        void set_no_element () { m_element = std::nullopt; }

        bool forward () const { return m_direction == TrackingDirection::Forward; }
    };

} // namespace impactx

#endif // IMPACTX_TRACKING_STATE_H
