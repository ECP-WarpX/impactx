/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Edoardo Zoni, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_ACCESSORS_H
#define IMPACTX_ELEMENTS_MIXIN_ACCESSORS_H

#include "traits.H"

#include "elements/Handle.H"

#include <AMReX_REAL.H>

#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>


namespace impactx::elements
{
    /** @file
     *
     * Property access on the ``KnownElements`` variant.
     *
     * These free functions wrap the ``std::visit`` dispatch that would otherwise be spelled
     * out at every call site that only wants to read or write one element property. Getters
     * and setters are overloaded on arity and named after the property itself, mirroring the
     * element member functions: ``ds(element)`` reads, ``ds(element, 0.1)`` writes.
     *
     * Since ``KnownElements`` is a ``std::variant`` over types in this namespace, ADL finds
     * these unqualified; the explicit ``elements::`` spelling is preferred for clarity.
     *
     * Getters are total: every lattice element provides ``ds()``, ``nslice()``, ``name()``,
     * ``has_name()`` and a static ``type``. Setters are partial and throw ``std::runtime_error``
     * on elements that do not support them, @see mixin::has_settable_ds_v,
     * mixin::has_settable_nslice_v and mixin::is_named_v for the exact sets.
     *
     * @note This header includes ``elements/All.H`` and therefore sits *above* the element
     *       definitions. Do not include it from an element header or from another header in
     *       ``mixin/``. Those need the type-level traits instead, @see traits.H.
     */

    /** Return the segment length of an element
     *
     * @param element_variant the element to query
     * @return value in meters, zero for zero-length elements
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline amrex::ParticleReal
    ds (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) { return element.ds(); }, element_variant);
    }

    /** Set the segment length of an element
     *
     * @param[in,out] element_variant the element to modify
     * @param ds segment length in meters
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline void
    ds (T_Element & element_variant, amrex::ParticleReal ds)
    {
        elements::visit([ds](auto & element) {
            if constexpr (mixin::has_settable_ds_v<decltype(element)>) {
                element.m_ds = ds;
            } else {
                throw std::runtime_error(
                    "elements::ds: this element has no settable segment length.");
            }
        }, element_variant);
    }

    /** Return the number of slices used for the application of collective effects
     *
     * @param element_variant the element to query
     * @return positive integer, one for elements that cannot be sliced
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline int
    nslice (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) { return element.nslice(); }, element_variant);
    }

    /** Set the number of slices used for the application of collective effects
     *
     * This subdivides the element's transport: the element then applies ``ds/nslice`` per
     * push. @see ScopedNslice to do this temporarily.
     *
     * @param[in,out] element_variant the element to modify
     * @param nslice number of slices, a positive integer
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline void
    nslice (T_Element & element_variant, int nslice)
    {
        elements::visit([nslice](auto & element) {
            if constexpr (mixin::has_settable_nslice_v<decltype(element)>) {
                if (nslice <= 0) {
                    throw std::invalid_argument(
                        "elements::nslice: nslice must be > 0.");
                }
                element.m_nslice = nslice;
            } else {
                throw std::runtime_error(
                    "elements::nslice: this element cannot be sliced.");
            }
        }, element_variant);
    }

    /** Return the length of a single slice of an element
     *
     * @param element_variant the element to query
     * @return ``ds/nslice`` in meters, zero for zero-length elements
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline amrex::ParticleReal
    slice_ds (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) {
            return element.ds() / element.nslice();
        }, element_variant);
    }

    /** Return the name of an element
     *
     * @param element_variant the element to query
     * @return the element name, or an empty string if it has none
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline std::string
    name (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) {
            return element.has_name() ? element.name() : std::string{};
        }, element_variant);
    }

    /** Set the name of an element
     *
     * @param[in,out] element_variant the element to modify
     * @param name the new element name
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline void
    name (T_Element & element_variant, std::string const & name)
    {
        elements::visit([&name](auto & element) {
            if constexpr (mixin::is_named_v<decltype(element)>) {
                element.set_name(name);
            } else {
                throw std::runtime_error(
                    "elements::name: this element has no settable name.");
            }
        }, element_variant);
    }

    /** Does this element carry a name?
     *
     * @param element_variant the element to query
     * @return true if the element has a name
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline bool
    has_name (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) { return element.has_name(); }, element_variant);
    }

    /** Return the element type, e.g. "quad"
     *
     * @param element_variant the element to query
     * @return the element's static type name
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline std::string
    type (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) {
            return std::string(std::decay_t<decltype(element)>::type);
        }, element_variant);
    }

    /** Does this element have a finite length? (@see mixin::Thick)
     *
     * @param element_variant the element to query
     * @return true for finite-length elements
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline bool
    is_thick (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) {
            return mixin::is_thick_v<decltype(element)>;
        }, element_variant);
    }

    /** Does this element have zero length? (@see mixin::Thin)
     *
     * @param element_variant the element to query
     * @return true for zero-length elements
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline bool
    is_thin (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) {
            return mixin::is_thin_v<decltype(element)>;
        }, element_variant);
    }

    /** Can this element be subdivided by overriding its slice count?
     *
     * True for the elements that derive their per-push transport length from their own
     * ``nslice()``, so that raising the slice count shortens one push. Note that this is
     * narrower than "has a slice count": @see mixin::can_slice_v for which elements are
     * excluded and why.
     *
     * @param element_variant the element to query
     * @return true if the element can be subdivided, @see ScopedNslice
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline bool
    can_slice (T_Element const & element_variant)
    {
        return elements::visit([](auto const & element) {
            return mixin::can_slice_v<decltype(element)>;
        }, element_variant);
    }

    /** Reverse an element in place
     *
     * After this call, pushing particles through the element reverses the effect of the
     * original element. This is an involution.
     *
     * @param[in,out] element_variant the element to reverse
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline void
    reverse (T_Element & element_variant)
    {
        elements::visit([](auto & element) { element.reverse(); }, element_variant);
    }

    /** Close and deallocate all data and handles of an element
     *
     * @param[in,out] element_variant the element to finalize
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    inline void
    finalize (T_Element & element_variant)
    {
        elements::visit([](auto & element) { element.finalize(); }, element_variant);
    }

    /** Temporarily override an element's slice count, restoring it on scope exit
     *
     * Subdividing an element is how a caller shortens its transport: with the slice count
     * temporarily doubled, each push advances ``ds/(2*nslice)`` instead of ``ds/nslice``.
     * Restoring on scope exit keeps the lattice consistent even if the pushes in between
     * throw.
     *
     * The element must be sliceable, @see can_slice.
     */
    template<typename T_Element, std::enable_if_t<is_element_repr_v<T_Element>, int> = 0>
    class ScopedNslice
    {
    public:
        /** Override an element's slice count for the lifetime of this object
         *
         * @param[in,out] element_variant the element to subdivide
         * @param nslice the temporary number of slices, a positive integer
         */
        ScopedNslice (T_Element & element_variant, int nslice)
        : m_element_variant(element_variant),
          m_nslice_original(elements::nslice(element_variant))
        {
            if (!elements::can_slice(m_element_variant)) {
                throw std::runtime_error(
                    "elements::ScopedNslice: this element cannot be sliced.");
            }
            elements::nslice(m_element_variant, nslice);
        }

        ~ScopedNslice ()
        {
            elements::nslice(m_element_variant, m_nslice_original);
        }

        ScopedNslice (ScopedNslice const &) = delete;
        ScopedNslice & operator= (ScopedNslice const &) = delete;
        ScopedNslice (ScopedNslice &&) = delete;
        ScopedNslice & operator= (ScopedNslice &&) = delete;

    private:
        T_Element & m_element_variant; //! the subdivided element
        int m_nslice_original; //! slice count to restore on scope exit
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_ACCESSORS_H
