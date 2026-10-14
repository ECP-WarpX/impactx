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

#include "elements/helper/Handle.H"

#include <atomic>
#include <cstddef>
#include <cstdint>
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
     *
     * The edits are defined in Lattice.cpp: each instantiates vector operations over the
     * variant of every element type, which would otherwise be compiled again in every file
     * that includes this header.
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

        /** A copy is a new lattice holding the same elements
         *
         * The traversal count is deliberately not copied. It says something about a loop
         * walking *this* instance, not about which elements a lattice holds, so a copy
         * made while the original is being walked would refuse every structural edit for
         * the rest of its life.
         */
        Lattice (Lattice const & other);

        /** @see Lattice(Lattice const &); a lattice being walked cannot be replaced */
        Lattice & operator= (Lattice const & other);

        /** @see Lattice(Lattice const &) */
        Lattice (Lattice && other) noexcept;

        /** @see Lattice(Lattice const &); a lattice being walked cannot be replaced */
        Lattice & operator= (Lattice && other) noexcept(false);

        /** Take shared ownership of an element value
         *
         * An rvalue is moved into the lattice, an lvalue copied.
         *
         * @param element the element to store
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
        void push_back (elements::ElementHandle handle);

        /** Insert an element before a position
         *
         * @param i the position to insert before, at most @see size
         * @param handle a non-null element handle
         */
        void insert (size_type i, elements::ElementHandle handle);

        /** Replace the element at a position
         *
         * Other positions holding the element that was there are unaffected.
         *
         * @param i the position to replace, less than @see size
         * @param handle a non-null element handle
         */
        void replace (size_type i, elements::ElementHandle handle);

        /** Remove the element at a position
         *
         * Other positions holding the same element are unaffected: this removes one
         * occurrence, not the element.
         */
        void erase (size_type i);

        /** Remove the last element */
        void pop_back ();

        /** Remove all elements */
        void clear ();

        [[nodiscard]] size_type size () const noexcept { return m_elements.size(); }
        [[nodiscard]] bool empty () const noexcept { return m_elements.empty(); }

        [[nodiscard]] elements::ElementHandle & operator[] (size_type i) { return m_elements[i]; }
        [[nodiscard]] elements::ElementHandle const & operator[] (size_type i) const { return m_elements[i]; }

        [[nodiscard]] elements::ElementHandle & at (size_type i) { return m_elements.at(i); }
        [[nodiscard]] elements::ElementHandle const & at (size_type i) const { return m_elements.at(i); }

        [[nodiscard]] elements::ElementHandle & front ()
        {
            require_not_empty("read the first element of");
            return m_elements.front();
        }
        [[nodiscard]] elements::ElementHandle const & front () const
        {
            require_not_empty("read the first element of");
            return m_elements.front();
        }
        [[nodiscard]] elements::ElementHandle & back ()
        {
            require_not_empty("read the last element of");
            return m_elements.back();
        }
        [[nodiscard]] elements::ElementHandle const & back () const
        {
            require_not_empty("read the last element of");
            return m_elements.back();
        }

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
        void assign (storage_type elements);

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
        /** Throw unless @p i names an element that is there
         *
         * @param i the position that was asked for
         * @param what the operation that was attempted, used in the error message
         */
        void require_element_at (size_type i, char const * what) const;

        /** Throw unless @p i is a position an element can be inserted before
         *
         * One past the last element is a valid place to insert.
         *
         * @param i the position that was asked for
         */
        void require_insertion_point (size_type i) const;

        /** Throw if there is no element to operate on
         *
         * @param what the operation that was attempted, used in the error message
         */
        void require_not_empty (char const * what) const;

        /** Throw if a tracking traversal is walking this lattice
         *
         * @param what the edit that was attempted, used in the error message
         */
        void require_no_traversal (char const * what) const;

        storage_type m_elements; //! elements in lattice order; one may occupy several positions
        std::uint64_t m_generation = 0; //! count of structural edits, @see generation
        //! live @see TraversalGuard count; atomic, because an edit can come from another thread
        mutable std::atomic<std::size_t> m_active_traversals{0};
    };

} // namespace impactx

#endif // IMPACTX_LATTICE_H
