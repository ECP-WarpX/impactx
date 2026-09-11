/* Copyright 2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIAGNOSTICS_LINEAR_MAP_H
#define IMPACTX_DIAGNOSTICS_LINEAR_MAP_H

#include "elements/All.H"
#include "Lattice.H"
#include "particles/CovarianceMatrix.H"
#include "particles/ReferenceParticle.H"

#include <list>
#include <string>
#include <vector>


namespace impactx::diagnostics
{
    /** Policy for how @ref map_trace and @ref linear_map handle elements
     *  whose @c transport_map evaluation throws, i.e. elements that do not
     *  provide an analytic linear transport map (for example
     *  @c NonlinearLens or @c Programmable).
     *
     *  The distinction between the two identity-substituting policies is
     *  deliberate: @c IdentityWithWarning is appropriate when the user has
     *  not asked for a fallback (e.g. from the Twiss diagnostic) and should
     *  be made aware that part of the lattice was ignored;
     *  @c IdentitySilent is appropriate when the user has explicitly opted
     *  in to the fallback (e.g. through @c transfer_map's
     *  @c fallback_identity_map=True kwarg) and does not want repeated
     *  warnings.
     */
    enum class OnMissingLinearMap
    {
        IdentityWithWarning,  //!< substitute the identity map, record one warning per element type
        IdentitySilent,       //!< substitute the identity map without emitting a warning
        Throw                 //!< re-raise as @c std::runtime_error naming the offending element
    };

    /** One entry in a per-element trace of the cumulative linear transport map.
     *
     *  Each entry records the cumulative 6x6 linear transport map from the
     *  start of the traversal (where the map is the identity) up to the exit
     *  of one lattice element, along with a small amount of metadata
     *  describing that element. The entries in @c map_trace form a sequence
     *  parallel to the element list, with one extra leading entry for the
     *  starting position.
     *
     *  Phase-space ordering in @c M_cumulative is (x, px, y, py, t, pt), with
     *  x, y transverse positions (meters), px, py transverse canonical momenta
     *  normalized by the reference momentum, t arrival-time deviation relative
     *  to the reference (in meters of c*dt), and pt the momentum deviation
     *  relative to the reference (dimensionless).
     */
    struct MapTraceEntry
    {
        amrex::ParticleReal s;          //!< integrated path length along the reference orbit, in meters
        std::string element_name;       //!< user-supplied element name; empty if the element was not named
        std::string element_type;       //!< element type string (e.g. "Drift", "Quad", "Sbend")
        Map6x6 M_cumulative;            //!< cumulative 6x6 linear transport map from the start up to the exit of this element
    };

    /** Compose the linear transport map element-by-element along a lattice.
     *
     *  Walks the lattice once, advancing a local copy of the reference
     *  particle through each element's non-linear reference-particle push and
     *  multiplying the cumulative 6x6 linear transport map by each element's
     *  per-slice @c transport_map. The caller's reference particle is not
     *  modified.
     *
     *  Slice ordering inside each element matches the envelope tracker in
     *  @c src/tracking/envelope.cpp : for each of @c nslice sub-steps, first
     *  advance the reference particle by one slice, then evaluate
     *  @c transport_map at the post-advance reference state. This keeps both
     *  code paths consistent for elements whose linear map depends on the
     *  reference momentum through the element (for example radio-frequency
     *  cavities).
     *
     *  Effects that are part of the non-linear tracker but not of the linear
     *  map are deliberately skipped here: space charge, Coherent Synchrotron
     *  Radiation (CSR), and Incoherent Synchrotron Radiation (ISR). The
     *  result is a purely element-by-element linear-optics trace appropriate
     *  for Twiss / tune / chromaticity / dispersion diagnostics.
     *
     *  Elements whose @c transport_map throws (for example @c NonlinearLens
     *  and @c Programmable explicitly reject a linear-map evaluation) are
     *  handled according to @p on_missing: either replaced by the identity
     *  (with a single warning per offending element type recorded via the
     *  ABLASTR warn manager) or re-raised as a @c std::runtime_error naming
     *  the offending element.
     *
     * @param[in] lattice         list of lattice elements to traverse
     * @param[in] ref_part_init   reference particle at the starting s
     * @param[in] on_missing      policy for elements without a linear map;
     *                            default is identity + warning
     * @return                    vector of entries, size @c lattice.size() + 1 ;
     *                            the first entry is the identity map at the
     *                            starting s, each subsequent entry is the
     *                            cumulative map at the exit of one element.
     */
    std::vector<MapTraceEntry>
    map_trace (
        Lattice const & lattice,
        RefPart const & ref_part_init,
        OnMissingLinearMap on_missing = OnMissingLinearMap::IdentityWithWarning
    );

    /** Return the end-to-end linear transport map of a lattice.
     *
     *  Convenience wrapper around @ref map_trace that returns only the final
     *  cumulative map. For a ring, this is the one-turn map used as input to
     *  the Wolski eigendecomposition [Wolski, Phys. Rev. ST Accel. Beams 9,
     *  024001 (2006); doi:10.1103/PhysRevSTAB.9.024001] of the periodic
     *  Twiss, tune, and dispersion diagnostics in
     *  @c impactx.twiss_lattice.compute_twiss (see the Python module for
     *  the full bibliography, including the Courant-Snyder [Ann. Phys. 3, 1
     *  (1958)] and Mais-Ripken [DESY M-82-05, DESY 88-114] references that
     *  underlie the coupled parameterization).
     *
     * @param[in] lattice         list of lattice elements
     * @param[in] ref_part_init   reference particle at the starting s
     * @param[in] on_missing      policy for elements without a linear map;
     *                            default is identity + warning
     * @return                    end-to-end 6x6 linear transport map
     */
    Map6x6
    linear_map (
        Lattice const & lattice,
        RefPart const & ref_part_init,
        OnMissingLinearMap on_missing = OnMissingLinearMap::IdentityWithWarning
    );

} // namespace impactx::diagnostics

#endif // IMPACTX_DIAGNOSTICS_LINEAR_MAP_H
