/* Copyright 2022-2026 The ImpactX Community
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PYTHON_LATTICEOWNERS_H
#define IMPACTX_PYTHON_LATTICEOWNERS_H

#include "pyImpactX.H"
#include "Lattice.H"

#include <elements/Handle.H>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>


namespace impactx::python
{
    namespace py = pybind11;

    /** Take a share of the element a Python object wraps.
     *
     * The returned handle refers to that very object, not to a copy of it, which is what
     * lets a lattice entry and a Python variable denote the same element.
     *
     * The dispatch is generated from ``KnownElements``: adding an element type to that
     * variant extends this automatically, and a type that is not in it cannot be inserted.
     */
    template<typename... T_Elements>
    elements::ElementHandle
    handle_from_python (py::handle obj, std::variant<T_Elements...> *)
    {
        elements::ElementHandle handle;
        bool matched = false;

        auto const try_alternative = [&](auto * tag) {
            using element_type = std::remove_pointer_t<decltype(tag)>;
            if (matched) { return; }
            if (py::isinstance<element_type>(obj))
            {
                handle = obj.template cast<std::shared_ptr<element_type>>();
                matched = true;
            }
        };
        (try_alternative(static_cast<T_Elements *>(nullptr)), ...);

        if (!matched)
        {
            throw py::type_error(
                "Lattice: expected a lattice element, got '" +
                std::string(py::str(py::type::of(obj).attr("__name__"))) + "'");
        }
        return handle;
    }

    /** @see handle_from_python */
    inline elements::ElementHandle handle_from_python (py::handle obj)
    {
        return handle_from_python(obj, static_cast<elements::KnownElements *>(nullptr));
    }

    /** Reject use of a view whose simulation is gone
     *
     * A lattice reached through ``sim.lattice`` is a member of that simulation, and the
     * view refers to it without keeping the simulation alive -- a strong link there would
     * make a reference cycle and defer the simulation's destruction past AMReX shutdown. So
     * the view can outlive what it points at, and says so rather than reading freed memory.
     *
     * Every binding that reaches the C++ lattice goes through @see lattice_of.
     */
    inline void require_parent_alive (py::object const & self)
    {
        if (!py::hasattr(self, "_parent_sim")) { return; }  // a standalone lattice
        if (self.attr("_parent_sim")().is_none())
        {
            throw std::runtime_error(
                "This lattice belonged to a simulation that no longer exists. Keep a "
                "reference to the simulation for as long as you use its lattice.");
        }
    }

    /** The C++ lattice behind a bound lattice object, once it is safe to touch
     *
     * @param self the Python lattice object
     */
    inline Lattice & lattice_of (py::object const & self)
    {
        require_parent_alive(self);
        return self.cast<Lattice &>();
    }

    /** The Python objects a lattice keeps alive, one entry per lattice position.
     *
     * Sharing a C++ element is not by itself enough to give a user back the object they
     * put in: a Python subclass, its ``__dict__`` and its callbacks live in the wrapper,
     * and pybind11 will mint a fresh, plain wrapper if the original is collected. The
     * lattice therefore holds a reference to each wrapper for as long as it holds the
     * element.
     *
     * The list lives in the lattice's instance dictionary rather than in the C++ class, so
     * it is visible to Python's cyclic garbage collector: a lattice holding an element
     * that refers back to the simulation still collects.
     *
     * An entry is ``None`` when the element did not come from Python -- parsed from an
     * input file, say. It is filled in the first time that position is handed out.
     */
    class Owners
    {
    public:
        /** Access the owner list of a bound lattice, reconciling it if needed
         *
         * @param self the Python lattice object
         * @param lattice its C++ lattice
         */
        Owners (py::object const & self, Lattice const & lattice)
          : m_self(self), m_lattice(lattice)
        {
            require_parent_alive(self);

            py::object list = py::none();
            if (py::hasattr(self, attr_name()))
            {
                list = self.attr(attr_name());
            }

            // A change made from C++ -- input parsing, element insertion, a filtered edit --
            // leaves this list describing a lattice that is no longer there. The lattice
            // counts its structural edits, so compare that count rather than the length: a
            // replacement of the same size changes no length but every element.
            bool const stale =
                !py::hasattr(self, generation_attr()) ||
                py::cast<std::uint64_t>(self.attr(generation_attr())) != lattice.generation();

            if (list.is_none() || stale || py::len(list) != lattice.size())
            {
                // Rebuild, keeping the wrappers whose element is still at the same position
                // so their identity survives.
                py::list rebuilt;
                py::list const previous = list.is_none() ? py::list() : py::cast<py::list>(list);
                for (Lattice::size_type i = 0; i < lattice.size(); ++i)
                {
                    py::object keep = py::none();
                    if (i < py::len(previous))
                    {
                        py::object const candidate = previous[i];
                        if (!candidate.is_none() && owns(candidate, lattice[i]))
                        {
                            keep = candidate;
                        }
                    }
                    rebuilt.append(keep);
                }
                self.attr(attr_name()) = rebuilt;
                list = rebuilt;
            }

            m_list = py::cast<py::list>(list);
        }

        /** Record which version of the lattice this owner list describes
         *
         * Anything this scope did to the lattice is done by now, so the count recorded here
         * is the one the list matches. If the two are out of step -- a mutation that threw
         * between changing the lattice and changing the list -- nothing is recorded, and the
         * next access rebuilds.
         */
        ~Owners ()
        {
            try
            {
                // Releasing a wrapper can run Python -- an element subclass' `__del__`, a
                // weakref callback -- which may reach the lattice again and install a list
                // of its own. Recording a version for a list that is no longer the lattice's
                // would make that state look current forever, so leave it unrecorded and let
                // the next access rebuild.
                bool const still_ours =
                    py::hasattr(m_self, attr_name()) && m_list.is(m_self.attr(attr_name()));

                if (still_ours && py::len(m_list) == m_lattice.size())
                {
                    m_self.attr(generation_attr()) = m_lattice.generation();
                }
            }
            catch (...) {}  // a destructor may not throw; a missing stamp only costs a rebuild
        }

        Owners (Owners const &) = delete;
        Owners & operator= (Owners const &) = delete;
        Owners (Owners &&) = delete;
        Owners & operator= (Owners &&) = delete;

        /** The exact Python object for a lattice position
         *
         * If the element did not come from Python, a wrapper is created once and kept, so
         * that repeated access returns the same object.
         */
        [[nodiscard]] py::object get (Lattice::size_type i)
        {
            py::object owner = m_list[i];
            if (owner.is_none())
            {
                // Hand out a share of the element itself. pybind11 returns the wrapper
                // already registered for this address if one exists, so this creates a new
                // object only the first time.
                owner = std::visit([](auto const & el) { return py::cast(el); }, m_lattice[i]);
                m_list[i] = owner;
            }
            return owner;
        }

        /** The recorded object for a position, without creating one
         *
         * ``None`` when the element did not come from Python and has not been handed out
         * yet. Use this to carry entries across a rebuild: @see get would materialize a
         * wrapper for every position it is asked about.
         */
        [[nodiscard]] py::object owner_at (Lattice::size_type i) const
        {
            return m_list[i];
        }

        /** Replace the whole list in one step
         *
         * Rebuilding is how the operations that change many positions at once stay linear:
         * removing or inserting one position at a time moves the tail of the list each
         * time.
         *
         * @param list one entry per position of the lattice as it now is
         */
        void assign (py::list list)
        {
            m_list = std::move(list);
            m_self.attr(attr_name()) = m_list;
        }

        /** Record the object owning the element just appended */
        void append (py::handle obj) { m_list.append(obj); }

        /** Record the object owning the element at a position */
        void set (Lattice::size_type i, py::handle obj) { m_list[i] = obj; }

        /** Record the object owning an element inserted at a position */
        void insert (Lattice::size_type i, py::handle obj)
        {
            m_list.attr("insert")(i, obj);
        }

        /** Forget the owner of a position */
        void erase (Lattice::size_type i)
        {
            m_list.attr("pop")(i);
        }

        /** Forget every owner */
        void clear ()
        {
            m_list = py::list();
            m_self.attr(attr_name()) = m_list;
        }

        [[nodiscard]] py::list const & list () const noexcept { return m_list; }

    private:
        /** Does this Python object wrap the element this handle owns? */
        static bool owns (py::handle obj, elements::ElementHandle const & handle)
        {
            try
            {
                return elements::address(handle_from_python(obj)) == elements::address(handle);
            }
            catch (py::type_error const &)
            {
                return false;
            }
        }

        static char const * attr_name () { return "_element_owners"; }
        static char const * generation_attr () { return "_element_owners_generation"; }

        py::object m_self;
        Lattice const & m_lattice;
        py::list m_list;
    };

} // namespace impactx::python

#endif // IMPACTX_PYTHON_LATTICEOWNERS_H
