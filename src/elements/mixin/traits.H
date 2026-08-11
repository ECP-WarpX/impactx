/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Edoardo Zoni, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_TRAITS_H
#define IMPACTX_ELEMENTS_MIXIN_TRAITS_H

#include "alignment.H"
#include "named.H"
#include "pipeaperture.H"
#include "spintransport.H"
#include "thick.H"
#include "thin.H"

#include <AMReX_REAL.H>

#include <type_traits>
#include <utility>


namespace impactx::elements::mixin
{
    /** @file
     *
     * Compile-time traits over a *concrete* lattice element type.
     *
     * These replace the hand-spelled
     * ``std::is_base_of_v<mixin::X, std::remove_cv_t<std::remove_reference_t<decltype(el)>>>``
     * at the call sites: the decay is applied inside the trait, so a trait can be used
     * directly on the deduced type of a ``std::visit`` or ``BeamOptic`` lambda parameter.
     *
     * This header deliberately does not include ``elements/All.H``, so it can be used from
     * element headers, from other mixins (e.g. @see beamoptic.H) and from the Python
     * bindings. For dispatch over the ``KnownElements`` variant @see accessors.H.
     */

    /** Does this element have a finite length (@see Thick)? */
    template<typename T_Element>
    inline constexpr bool is_thick_v = std::is_base_of_v<Thick, std::decay_t<T_Element>>;

    /** Does this element have zero length (@see Thin)? */
    template<typename T_Element>
    inline constexpr bool is_thin_v = std::is_base_of_v<Thin, std::decay_t<T_Element>>;

    /** Does this element carry a user-settable name (@see Named)? */
    template<typename T_Element>
    inline constexpr bool is_named_v = std::is_base_of_v<Named, std::decay_t<T_Element>>;

    /** Does this element support alignment errors (@see Alignment)? */
    template<typename T_Element>
    inline constexpr bool has_alignment_v = std::is_base_of_v<Alignment, std::decay_t<T_Element>>;

    /** Does this element support a transverse beam pipe aperture (@see PipeAperture)? */
    template<typename T_Element>
    inline constexpr bool has_pipe_aperture_v = std::is_base_of_v<PipeAperture, std::decay_t<T_Element>>;

    /** Does this element support spin transport (@see SpinTransport)? */
    template<typename T_Element>
    inline constexpr bool has_spin_transport_v = std::is_base_of_v<SpinTransport, std::decay_t<T_Element>>;

    namespace detail
    {
        /** Does this element have a publicly writable ``m_ds`` member?
         *
         * std::void_t detection idiom. Note that presence alone does not imply that ``m_ds``
         * is the element's segment length, @see has_settable_ds_v.
         */
        template<typename T_Element, typename = void>
        struct has_m_ds : std::false_type {};

        template<typename T_Element>
        struct has_m_ds<
            T_Element,
            std::void_t<decltype(std::declval<std::decay_t<T_Element> &>().m_ds = amrex::ParticleReal(0))>
        > : std::true_type {};

        /** Does this element have a publicly writable ``m_nslice`` member? */
        template<typename T_Element, typename = void>
        struct has_m_nslice : std::false_type {};

        template<typename T_Element>
        struct has_m_nslice<
            T_Element,
            std::void_t<decltype(std::declval<std::decay_t<T_Element> &>().m_nslice = 1)>
        > : std::true_type {};

    } // namespace detail

    /** Can this element's segment length be set?
     *
     * True for all finite-length elements: everything deriving from @see Thick, plus
     * Programmable, LinearMap and SpinMap, which carry their own ``m_ds``.
     *
     * The ``!is_thin_v`` term is load-bearing and *not* redundant with the member detection:
     * ThinDipole is a zero-length element (@see Thin, so ``ds()`` returns 0) that nevertheless
     * owns an ``m_ds`` member holding its *effective (equivalent) arc length* ``theta * rc``.
     * Detection alone would let a setter silently overwrite that with a segment length.
     */
    template<typename T_Element>
    inline constexpr bool has_settable_ds_v = std::conjunction_v<
        std::bool_constant<!is_thin_v<T_Element>>,
        detail::has_m_ds<T_Element>
    >;

    /** Can this element's slice count be set?
     *
     * True for everything deriving from @see Thick, plus Programmable, which carries its own
     * ``m_nslice``. False for LinearMap and SpinMap: those have a finite ``ds`` but their
     * ``nslice()`` returns a literal 1.
     */
    template<typename T_Element>
    inline constexpr bool has_settable_nslice_v = std::conjunction_v<
        std::bool_constant<!is_thin_v<T_Element>>,
        detail::has_m_nslice<T_Element>
    >;

    /** Can this element be subdivided by overriding its slice count?
     *
     * True exactly for the elements that derive their per-push transport length from their
     * own ``nslice()`` at push time, i.e. everything deriving from @see Thick. Only for those
     * does temporarily raising the slice count actually shorten the transport of one push,
     * @see ScopedNslice in helper/Accessors.H.
     *
     * This is deliberately narrower than @see has_settable_nslice_v:
     * - Programmable *stores* a slice count, but its push is an opaque user callback. The
     *   callback reads ``ds`` and ``nslice`` from the user's own Python object, while the
     *   lattice holds a copy — a temporary override is invisible to it, and the callback
     *   would transport a full slice on each of the two half-maps.
     * - LinearMap and SpinMap have a finite ``ds`` but no slice count to raise: their
     *   ``nslice()`` returns a literal 1.
     */
    template<typename T_Element>
    inline constexpr bool can_slice_v = is_thick_v<T_Element>;

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_TRAITS_H
