/* Copyright 2022-2026 The ImpactX Community
 *
 * Authors: Axel Huebl, Eric G. Stern
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"
#include "Lattice.H"
#include "LatticeOwners.H"

#include <diagnostics/LinearMap.H>
#include <elements/All.H>
#include <elements/mixin/lineartransport.H>
#include <elements/transformation/Insert.H>
#include <particles/CovarianceMatrix.H>
#include <particles/ReferenceParticle.H>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace py = pybind11;
using namespace impactx;

// at namespace scope: MSVC does not find a block-scope using-declaration of a function
// from inside a nested lambda, and every use below sits in one
using impactx::python::handle_from_python;
using impactx::python::lattice_of;


void init_lattice(py::module& me)
{
    using elements::KnownElements;
    namespace ix_diag = ::impactx::diagnostics;

    // The lattice holds the element objects it is given, as a Python list holds its items:
    // one object may sit at several positions, and changing it changes what is tracked.
    using KnownElementsList = Lattice;
    using impactx::python::Owners;

    // dynamic_attr: the lattice keeps the Python wrapper of each element alive in its
    // instance dictionary, where the cyclic garbage collector can see it.
    py::class_<KnownElementsList, py::smart_holder> kel(me, "KnownElementsList", py::dynamic_attr());

    /** The element a Python object refers to, as an identity key
     *
     * ``nullptr`` when the object is not a lattice element at all, so that asking whether
     * a string is in the lattice answers no rather than raising.
     *
     * Comparing elements rather than the Python objects wrapping them avoids creating a
     * wrapper for every position looked at: a search that finds nothing would otherwise
     * leave one behind for the whole lattice.
     */
    auto const element_address = [](py::handle obj) -> void const * {
        try
        {
            return elements::address(handle_from_python(obj));
        }
        catch (py::type_error const &)
        {
            return nullptr;
        }
    };

    /** Normalize and bounds-check a Python index, allowing negative indexing */
    auto const checked_index = [](KnownElementsList const & v, py::ssize_t index) -> KnownElementsList::size_type {
        py::ssize_t const n = static_cast<py::ssize_t>(v.size());
        if (index < 0) { index += n; }
        if (index < 0 || index >= n) { throw py::index_error("Index out of range"); }
        return static_cast<KnownElementsList::size_type>(index);
    };

    kel
        .def(py::init<>())

        .def("append",
             [](py::object self, py::object el) {
                 auto & v = lattice_of(self);
                 auto handle = handle_from_python(el);
                 Owners owners(self, v);
                 v.push_back(std::move(handle));
                 owners.append(el);
             },
             py::arg("element"),
             "Add an element to the end of the lattice.\n\n"
             "The lattice holds this element, so changing it afterwards changes what is\n"
             "tracked. Adding the same element twice places it at two positions. For two\n"
             "independent elements, add ``element.copy()`` or construct a second one."
        )

        .def("extend",
             [](py::object self, py::iterable const & l) {
                 auto & v = lattice_of(self);

                 // convert everything first, so a bad entry leaves the lattice unchanged
                 std::vector<std::pair<elements::ElementHandle, py::object>> pending;
                 for (auto const & item : l)
                 {
                     py::object el = py::reinterpret_borrow<py::object>(item);
                     pending.emplace_back(handle_from_python(el), el);
                 }

                 Owners owners(self, v);
                 for (auto & [handle, el] : pending)
                 {
                     v.push_back(std::move(handle));
                     owners.append(el);
                 }
                 return self;
             },
             py::arg("elements"),
             "Add several elements to the lattice."
        )

        .def_property_readonly("generation",
             [](py::object self) { return lattice_of(self).generation(); },
             "How often the sequence of elements changed.\n\n"
             "Counts structural edits only: changing a parameter on an element that is\n"
             "already in the lattice does not move anything and does not change this.\n"
             "A view that remembers positions compares this to notice they went stale.\n\n"
             "Compare it for equality, not by how far it moved: one call can edit the\n"
             "sequence in several steps, as deleting a selection does."
        )

        .def("size", [](py::object self) { return lattice_of(self).size(); })
        .def("is_empty", [](py::object self) { return lattice_of(self).empty(); })
        .def("__len__", [](py::object self) { return lattice_of(self).size(); },
             "The length of the list.")

        .def("clear",
             [](py::object self) {
                 auto & v = lattice_of(self);
                 Owners owners(self, v);
                 v.clear();
                 owners.clear();
             },
             "Remove all elements from the lattice."
        )

        .def("pop_back",
             [](py::object self) {
                 auto & v = lattice_of(self);
                 if (v.empty()) { throw py::index_error("pop from empty lattice"); }
                 Owners owners(self, v);
                 py::object const last = owners.get(v.size() - 1);
                 v.pop_back();
                 owners.erase(v.size());
                 return last;
             },
             "Remove and return the last element of the lattice."
        )

        .def("__iter__",
             [](py::object self) {
                 auto & v = lattice_of(self);
                 Owners owners(self, v);
                 // materialize the exact wrappers, so iteration yields the objects that
                 // were put in rather than fresh views of them
                 py::list out;
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i) { out.append(owners.get(i)); }
                 return out.attr("__iter__")();
             },
             "Iterate over the elements of the lattice."
        )

        .def("__getitem__",
             [checked_index](py::object self, py::ssize_t index) {
                 auto & v = lattice_of(self);
                 auto const i = checked_index(v, index);
                 Owners owners(self, v);
                 return owners.get(i);
             },
             py::arg("index"),
             "Return the element at a position."
        )

        .def("__setitem__",
             [checked_index](py::object self, py::ssize_t index, py::object el) {
                 auto & v = lattice_of(self);
                 auto const i = checked_index(v, index);
                 auto handle = handle_from_python(el);
                 Owners owners(self, v);
                 v.replace(i, std::move(handle));
                 owners.set(i, el);
             },
             py::arg("index"), py::arg("element"),
             "Replace the element at a position."
        )

        .def("__getitem__",
             [](py::object self, py::slice const & slice) {
                 auto & v = lattice_of(self);
                 size_t start = 0, stop = 0, step = 0, length = 0;
                 if (!slice.compute(v.size(), &start, &stop, &step, &length)) {
                     throw py::error_already_set();
                 }
                 Owners owners(self, v);

                 // A slice is a new lattice over the same elements, like a Python list
                 // slice. It is a plain lattice even when sliced from a subclass, as
                 // slicing a `list` subclass gives a plain `list`: the subclass may take
                 // constructor arguments we have nothing to pass, and its `extend` may
                 // mean something of its own.
                 py::object out = py::type::of<KnownElementsList>()();
                 auto & sliced = out.cast<KnownElementsList &>();
                 Owners sliced_owners(out, sliced);
                 for (size_t i = 0; i < length; ++i)
                 {
                     auto const at = start + i * step;
                     sliced.push_back(v[at]);
                     sliced_owners.append(owners.get(at));
                 }
                 return out;
             },
             py::arg("slice"),
             "Return a new lattice holding the selected elements."
        )

        .def("__setitem__",
             [](py::object self, py::slice const & slice, py::iterable const & value) {
                 auto & v = lattice_of(self);

                 // Convert first, for two reasons: a bad entry must leave the lattice
                 // untouched, and consuming the iterable can itself change the lattice --
                 // a generator is free to append to it. Positions computed before that
                 // would name the wrong elements, so the bounds are taken afterwards,
                 // against the length the lattice actually has. This is also the order
                 // list.__setitem__ uses.
                 std::vector<std::pair<elements::ElementHandle, py::object>> pending;
                 for (auto const & item : value)
                 {
                     py::object el = py::reinterpret_borrow<py::object>(item);
                     pending.emplace_back(handle_from_python(el), el);
                 }

                 size_t start = 0, stop = 0, step = 0, length = 0;
                 if (!slice.compute(v.size(), &start, &stop, &step, &length)) {
                     throw py::error_already_set();
                 }

                 if (step != 1 && pending.size() != length)
                 {
                     throw py::value_error(
                         "attempt to assign sequence of size " + std::to_string(pending.size()) +
                         " to extended slice of size " + std::to_string(length));
                 }

                 // replacing nothing with nothing changes nothing
                 if (length == 0 && pending.empty()) { return; }

                 // Replacing a position lets go of what was there: the C++ handle, and
                 // with the owner list the last reference to its Python wrapper. Either
                 // can run a subclass' `__del__` from inside the update, and a finalizer
                 // that reads the lattice sees a generation that has already moved, so it
                 // rebuilds the owner list -- leaving the writes below to land in a list
                 // that is no longer the lattice's. Hold both the displaced handles and
                 // their wrappers until every position and every owner is committed.
                 // Declared before `owners` so they outlive it and are released only once
                 // the owner list has been stamped.
                 std::vector<elements::ElementHandle> displaced;
                 py::list displaced_owners;
                 displaced.reserve(length);
                 for (size_t i = 0; i < length; ++i)
                 {
                     displaced.push_back(v[start + i * step]);
                 }

                 Owners owners(self, v);
                 for (size_t i = 0; i < length; ++i)
                 {
                     displaced_owners.append(owners.owner_at(start + i * step));
                 }
                 if (step == 1)
                 {
                     // A contiguous slice may change the length, so build the result in one
                     // pass: head, the new elements, tail. Erasing and inserting one position
                     // at a time moves everything after `start` on every step, which is
                     // quadratic -- and `lattice[:] = new` is the ordinary way to say
                     // "replace the lattice".
                     Lattice::storage_type next;
                     next.reserve(v.size() - length + pending.size());
                     py::list next_owners;

                     for (Lattice::size_type i = 0; i < start; ++i)
                     {
                         next.push_back(v[i]);
                         next_owners.append(owners.owner_at(i));
                     }
                     for (auto & [handle, el] : pending)
                     {
                         next.push_back(std::move(handle));
                         next_owners.append(el);
                     }
                     for (Lattice::size_type i = start + length; i < v.size(); ++i)
                     {
                         next.push_back(v[i]);
                         next_owners.append(owners.owner_at(i));
                     }

                     v.assign(std::move(next));
                     owners.assign(next_owners);
                 }
                 else
                 {
                     for (size_t i = 0; i < length; ++i)
                     {
                         auto const at = start + i * step;
                         v.replace(at, std::move(pending[i].first));
                         owners.set(at, pending[i].second);
                     }
                 }
             },
             py::arg("slice"), py::arg("elements"),
             "Replace the selected elements."
        )

        .def("__delitem__",
             [](py::object self, py::slice const & slice) {
                 auto & v = lattice_of(self);
                 size_t start = 0, stop = 0, step = 0, length = 0;
                 if (!slice.compute(v.size(), &start, &stop, &step, &length)) {
                     throw py::error_already_set();
                 }
                 // Nothing selected changes nothing, and must not count as an edit:
                 // that would void every selection taken on this lattice.
                 if (length == 0) { return; }

                 // Mark what goes, then keep the rest in one pass. Removing positions one
                 // at a time moves the tail of both the lattice and the owner list on every
                 // removal, which is quadratic in the length of the lattice; a slice over a
                 // long beamline is exactly where that bites. `step` is unsigned, so a
                 // negative step counts down through wraparound -- marking sidesteps the
                 // question of which position is the largest.
                 std::vector<bool> dropped(v.size(), false);
                 for (size_t i = 0; i < length; ++i)
                 {
                     dropped[start + i * step] = true;
                 }

                 Owners owners(self, v);

                 Lattice::storage_type kept;
                 kept.reserve(v.size() - length);
                 py::list kept_owners;
                 for (Lattice::size_type i = 0; i < v.size(); ++i)
                 {
                     if (dropped[i]) { continue; }
                     kept.push_back(v[i]);
                     kept_owners.append(owners.owner_at(i));
                 }

                 v.assign(std::move(kept));
                 owners.assign(kept_owners);
             },
             py::arg("slice"),
             "Remove the selected elements."
        )

        .def("__delitem__",
             [checked_index](py::object self, py::ssize_t index) {
                 auto & v = lattice_of(self);
                 auto const i = checked_index(v, index);
                 Owners owners(self, v);
                 v.erase(i);
                 owners.erase(i);
             },
             py::arg("index"),
             "Remove the element at a position.\n\n"
             "An element that also sits at other positions keeps those."
        )

        .def("insert",
             [](py::object self, py::ssize_t index, py::object el) {
                 auto & v = lattice_of(self);
                 auto handle = handle_from_python(el);

                 // Python list semantics: an out-of-range index clamps rather than raises
                 py::ssize_t const n = static_cast<py::ssize_t>(v.size());
                 if (index < 0) { index += n; }
                 index = std::clamp<py::ssize_t>(index, 0, n);

                 Owners owners(self, v);
                 v.insert(static_cast<KnownElementsList::size_type>(index), std::move(handle));
                 owners.insert(static_cast<KnownElementsList::size_type>(index), el);
             },
             py::arg("index"), py::arg("element"),
             "Insert an element before a position."
        )

        // By identity, deliberately, and by identity only. Elements compare by value, so
        // two distinct elements with the same parameters are equal; for locating an
        // occurrence in a lattice the caller means the element they hold, not another one
        // configured the same way. This is why ``element in lattice`` can be False while
        // ``lattice == [element]`` is True.
        .def("index",
             [element_address](py::object self, py::object el) {
                 auto & v = lattice_of(self);
                 auto const * wanted = element_address(el);
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                 {
                     if (wanted != nullptr && elements::address(v[i]) == wanted) { return i; }
                 }
                 throw py::value_error("element is not in the lattice");
             },
             py::arg("element"),
             "Return the first position holding this element."
        )

        .def("count",
             [element_address](py::object self, py::object el) {
                 auto & v = lattice_of(self);
                 auto const * wanted = element_address(el);
                 KnownElementsList::size_type n = 0;
                 if (wanted == nullptr) { return n; }
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                 {
                     if (elements::address(v[i]) == wanted) { ++n; }
                 }
                 return n;
             },
             py::arg("element"),
             "Return how many positions this element occupies."
        )

        .def("__contains__",
             [element_address](py::object self, py::object el) {
                 auto & v = lattice_of(self);
                 auto const * wanted = element_address(el);
                 if (wanted == nullptr) { return false; }
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                 {
                     if (elements::address(v[i]) == wanted) { return true; }
                 }
                 return false;
             },
             py::arg("element")
        )

        .def("remove",
             [element_address](py::object self, py::object el) {
                 auto & v = lattice_of(self);
                 auto const * wanted = element_address(el);
                 if (wanted != nullptr)
                 {
                     for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                     {
                         if (elements::address(v[i]) != wanted) { continue; }
                         Owners owners(self, v);
                         v.erase(i);
                         owners.erase(i);
                         return;
                     }
                 }
                 throw py::value_error("element is not in the lattice");
             },
             py::arg("element"),
             "Remove the first occurrence of this element."
        )

        .def("__reversed__",
             [](py::object self) {
                 auto & v = lattice_of(self);
                 Owners owners(self, v);
                 py::list out;
                 for (KnownElementsList::size_type i = v.size(); i-- > 0; ) { out.append(owners.get(i)); }
                 return out.attr("__iter__")();
             },
             "Iterate over the elements from the end."
        )

        .def(
            "transfer_map",
            // The reference particle is taken by value. pybind11 copies it
            // from the caller on entry, so the caller's ``sim.beam.ref`` is
            // not modified when the internal traversal advances the
            // reference.
            [](
                py::object self,
                RefPart ref,
                std::string order,
                bool fallback_identity_map
            )
            {
                auto const & v = lattice_of(self);
                if (order != "linear") {
                    throw std::runtime_error(
                        "Only the calculation of linear transfer maps is "
                        "currently supported."
                    );
                }
                auto const on_missing = fallback_identity_map
                    ? ix_diag::OnMissingLinearMap::IdentitySilent
                    : ix_diag::OnMissingLinearMap::Throw;
                return ix_diag::linear_map(v, ref, on_missing);
            },
            py::arg("ref"),
            py::arg("order") = "linear",
            py::arg("fallback_identity_map") = false,
            "Calculate the end-to-end transfer map of the elements in the list.\n\n"
            "Currently only the linear transfer map is implemented (``order=\"linear\"``);\n"
            "the ``order`` parameter is reserved for future higher-order extensions.\n"
            "In linear mode the 6x6 map is composed element by element, using each\n"
            "element's analytic per-slice linear transport map.\n\n"
            "Collective effects like space charge, Coherent/Incoherent Synchrotron\n"
            "Radiation (CSR/ISR), and wakefield effects are not applied here; the\n"
            "returned map describes the purely linear single-particle dynamics of the\n"
            "design lattice.\n\n"
            "Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).\n\n"
            ":param ref: reference particle at the starting s\n"
            ":param order: So far, only the calculation of linear transfer maps is supported.\n"
            ":param fallback_identity_map: For elements with an undefined transfer map in the lattice, assume the identity matrix.\n"
        )

        .def(
            "map_trace",
            [](py::object self, RefPart ref)  // intentional copy of ref
            {
                auto const & v = lattice_of(self);
                auto const trace = ix_diag::map_trace(v, ref);

                py::list out;
                for (auto const & e : trace)
                {
                    py::dict d;
                    d["s"] = e.s;
                    d["name"] = e.element_name;
                    d["type"] = e.element_type;
                    d["M"] = e.M_cumulative;
                    out.append(std::move(d));
                }
                return out;
            },
            py::arg("ref"),
            "Trace the cumulative 6x6 linear transport map element by element.\n\n"
            "The reference particle is passed by value (intentional copy); the\n"
            "caller's reference particle is not modified in place. This matches\n"
            "the convention used by ``transfer_map``.\n\n"
            "This per-element trace is what ``sim.twiss()`` consumes to transport\n"
            "Twiss functions through the lattice.\n\n"
            "If you only need the final cumulative map at the lattice exit,\n"
            "prefer ``transfer_map(ref)`` instead of indexing the last entry\n"
            "of ``map_trace(ref)``.\n\n"
            ":param ref: A reference particle.\n"
            ":return: A list of dictionaries, one per lattice element plus a\n"
            "    leading entry for the starting position. Each entry contains:\n\n"
            "    * ``s``    -- integrated path length along the reference orbit,\n"
            "      in meters;\n"
            "    * ``name`` -- user-supplied element name (empty string if not\n"
            "      named);\n"
            "    * ``type`` -- element type string (e.g. ``\"Drift\"``,\n"
            "      ``\"Quad\"``, ``\"Sbend\"``);\n"
            "    * ``M``    -- cumulative 6x6 linear transport map from the\n"
            "      start of the lattice to the exit of this element (a\n"
            "      ``Map6x6`` instance; call ``.to_numpy()`` for a standard\n"
            "      C-ordered NumPy array).\n\n"
            "    The first entry always has the identity map at the starting\n"
            "    ``s``; the last entry contains the same map as ``transfer_map``."
        )
    ;


    // lattice transformations
    py::module_ met = me.def_submodule(
        "transformation",
        "Transform and modify lattices"
    );

    met.def(
        "insert_element_every_ds",
        [](py::object source, amrex::ParticleReal ds, elements::KnownElements element) {
            auto & from = lattice_of(source);

            auto result = impactx::elements::transformation::insert_element_every_ds(
                from, ds, std::move(element));

            // Carry over the Python object of every element that came through unsplit, so
            // that it is still the element the caller put in -- with its subclass and its
            // attributes -- rather than a plain wrapper minted on first access. Positions
            // holding a split half, or a copy of the inserted element, have no object of
            // their own and are filled in when they are first asked for.
            py::object out = py::type::of<KnownElementsList>()();
            auto & into = out.cast<KnownElementsList &>();
            Owners from_owners(source, from);
            Owners into_owners(out, into);

            // Where each element of the source sits, so that a position holding a newly
            // made element -- an inserted copy, or half of one that was split -- does not
            // cost a scan and does not lose track of the elements after it. An element at
            // several positions is matched to them in order.
            std::unordered_map<void const *, std::vector<KnownElementsList::size_type>> positions;
            for (KnownElementsList::size_type i = 0; i < from.size(); ++i)
            {
                positions[elements::address(from[i])].push_back(i);
            }
            std::unordered_map<void const *, std::size_t> taken;

            for (KnownElementsList::size_type i = 0; i < result.size(); ++i)
            {
                into.push_back(result[i]);

                py::object owner = py::none();
                auto const * key = elements::address(result[i]);
                auto const found = positions.find(key);
                if (found != positions.end())
                {
                    auto & next = taken[key];
                    if (next < found->second.size())
                    {
                        owner = from_owners.owner_at(found->second[next]);
                        ++next;
                    }
                }
                into_owners.append(owner);
            }
            return out;
        },
        py::arg("list"),
        py::arg("ds"),
        py::arg("element"),
        "Insert an element every s into an element list.\n\n"
        "Returns a new lattice. Elements that are not split are the same objects as\n"
        "in the lattice given; an element that a split falls inside is replaced by two\n"
        "new elements covering its halves, which carry its parameters but not a Python\n"
        "subclass or attributes."
    );
}
