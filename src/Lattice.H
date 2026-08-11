/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_LATTICE_H
#define IMPACTX_LATTICE_H

#include "elements/Handle.H"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>


namespace impactx
{
    /** The sequence of elements a beam is tracked through.
     *
     * Entries are shared handles, so the same physical element can occupy several
     * positions, be shared between lattices, and outlive any one of them. Storage is
     * contiguous, which makes indexed access constant-time -- selection and filtering walk
     * the lattice by index.
     *
     * Structural edits validate their whole input before changing anything, so a rejected
     * edit leaves the lattice as it was rather than half-applied.
     */
    class Lattice
    {
    public:
        using storage_type = std::vector<elements::ElementHandle>;
        using value_type = storage_type::value_type;
        using size_type = storage_type::size_type;
        using iterator = storage_type::iterator;
        using const_iterator = storage_type::const_iterator;
        using reverse_iterator = storage_type::reverse_iterator;
        using const_reverse_iterator = storage_type::const_reverse_iterator;

        Lattice () = default;

        /** Take shared ownership of a freshly constructed element
         *
         * @tparam T_Element the concrete element type
         * @param args arguments forwarded to the element constructor
         * @return a handle to the element just added
         */
        template<typename T_Element, typename... T_Args>
        elements::ElementHandle & emplace (T_Args&&... args)
        {
            require_no_traversal("add an element");
            m_elements.push_back(elements::make_handle<T_Element>(std::forward<T_Args>(args)...));
            ++m_generation;
            return m_elements.back();
        }

        /** Take shared ownership of an element, copying the value given
         *
         * @param element the element to copy into the lattice
         */
        template<typename T_Element>
        void emplace_back (T_Element && element)
        {
            require_no_traversal("add an element");
            m_elements.push_back(
                elements::make_handle(elements::KnownElements(std::forward<T_Element>(element))));
            ++m_generation;
        }

        /** Share an element that is already owned elsewhere
         *
         * @param handle a non-null element handle
         */
        void push_back (elements::ElementHandle handle)
        {
            elements::require_valid(handle, "insert");
            require_no_traversal("add an element");
            m_elements.push_back(std::move(handle));
            ++m_generation;
        }

        /** Insert an element before a position
         *
         * @param i the position to insert before, at most @see size
         * @param handle a non-null element handle
         */
        void insert (size_type i, elements::ElementHandle handle)
        {
            elements::require_valid(handle, "insert");
            require_no_traversal("add an element");
            m_elements.insert(m_elements.begin() + static_cast<std::ptrdiff_t>(i), std::move(handle));
            ++m_generation;
        }

        /** Remove the element at a position
         *
         * Other positions holding the same element are unaffected: this removes one
         * occurrence, not the element.
         */
        void erase (size_type i)
        {
            require_no_traversal("remove an element");
            m_elements.erase(m_elements.begin() + static_cast<std::ptrdiff_t>(i));
            ++m_generation;
        }

        /** Remove the last element */
        void pop_back ()
        {
            require_no_traversal("remove an element");
            m_elements.pop_back();
            ++m_generation;
        }

        /** Remove all elements */
        void clear ()
        {
            require_no_traversal("clear the lattice");
            m_elements.clear();
            ++m_generation;
        }

        [[nodiscard]] size_type size () const noexcept { return m_elements.size(); }
        [[nodiscard]] bool empty () const noexcept { return m_elements.empty(); }

        [[nodiscard]] elements::ElementHandle & operator[] (size_type i) { return m_elements[i]; }
        [[nodiscard]] elements::ElementHandle const & operator[] (size_type i) const { return m_elements[i]; }

        [[nodiscard]] elements::ElementHandle & at (size_type i) { return m_elements.at(i); }
        [[nodiscard]] elements::ElementHandle const & at (size_type i) const { return m_elements.at(i); }

        [[nodiscard]] elements::ElementHandle & front () { return m_elements.front(); }
        [[nodiscard]] elements::ElementHandle const & front () const { return m_elements.front(); }
        [[nodiscard]] elements::ElementHandle & back () { return m_elements.back(); }
        [[nodiscard]] elements::ElementHandle const & back () const { return m_elements.back(); }

        [[nodiscard]] iterator begin () noexcept { return m_elements.begin(); }
        [[nodiscard]] iterator end () noexcept { return m_elements.end(); }
        [[nodiscard]] const_iterator begin () const noexcept { return m_elements.begin(); }
        [[nodiscard]] const_iterator end () const noexcept { return m_elements.end(); }
        [[nodiscard]] reverse_iterator rbegin () noexcept { return m_elements.rbegin(); }
        [[nodiscard]] reverse_iterator rend () noexcept { return m_elements.rend(); }
        [[nodiscard]] const_reverse_iterator rbegin () const noexcept { return m_elements.rbegin(); }
        [[nodiscard]] const_reverse_iterator rend () const noexcept { return m_elements.rend(); }

        /** Replace the whole sequence in one transaction
         *
         * Every handle is checked before anything is stored, so a null entry anywhere
         * leaves the lattice untouched.
         */
        void assign (storage_type elements)
        {
            for (auto const & handle : elements)
            {
                elements::require_valid(handle, "replace the lattice");
            }
            require_no_traversal("replace the lattice");
            m_elements = std::move(elements);
            ++m_generation;
        }

        /** How often the sequence of elements changed
         *
         * Views that remember positions compare this to detect that their indices went
         * stale. It counts structural edits only; changing a parameter on an element that
         * is already in the lattice does not move anything and does not bump it.
         */
        [[nodiscard]] std::uint64_t generation () const noexcept { return m_generation; }

        /** Marks a tracking traversal as in progress for its lifetime.
         *
         * While one exists, structural edits throw instead of invalidating the iteration
         * that is walking the lattice. Reads, and writes to an element's own parameters,
         * stay legal -- retuning an element from a tracking hook is a supported thing to do.
         */
        class TraversalGuard
        {
        public:
            explicit TraversalGuard (Lattice const & lattice)
              : m_lattice(&lattice)
            {
                ++m_lattice->m_active_traversals;
            }

            ~TraversalGuard () { --m_lattice->m_active_traversals; }

            TraversalGuard (TraversalGuard const &) = delete;
            TraversalGuard & operator= (TraversalGuard const &) = delete;
            TraversalGuard (TraversalGuard &&) = delete;
            TraversalGuard & operator= (TraversalGuard &&) = delete;

        private:
            Lattice const * m_lattice;
        };

        [[nodiscard]] bool is_being_traversed () const noexcept { return m_active_traversals > 0; }

    private:
        /** Throw if a tracking traversal is walking this lattice
         *
         * @param what the edit that was attempted, used in the error message
         */
        void require_no_traversal (char const * what) const
        {
            if (m_active_traversals > 0)
            {
                throw std::runtime_error(
                    std::string("Lattice: cannot ") + what + " while tracking through it. "
                    "Element parameters may be changed during tracking, but the sequence of "
                    "elements may not.");
            }
        }

        storage_type m_elements;
        std::uint64_t m_generation = 0;
        mutable std::size_t m_active_traversals = 0;
    };

} // namespace impactx

#endif // IMPACTX_LATTICE_H
