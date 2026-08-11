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
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace py = pybind11;
using namespace impactx;


void init_lattice(py::module& me)
{
    using elements::KnownElements;
    namespace ix_diag = ::impactx::diagnostics;

    // Elements are shared, not copied: appending keeps the object the caller passed in,
    // exactly as appending to a Python list does. Two positions holding the same object
    // are two occurrences of one element. Independent elements come from constructing one
    // per position, or from an explicit copy.
    using KnownElementsList = Lattice;
    using impactx::python::Owners;
    using impactx::python::handle_from_python;

    // dynamic_attr: the lattice keeps the Python wrapper of each element alive in its
    // instance dictionary, where the cyclic garbage collector can see it.
    py::class_<KnownElementsList, py::smart_holder> kel(me, "KnownElementsList", py::dynamic_attr());

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
                 auto & v = self.cast<KnownElementsList &>();
                 auto handle = handle_from_python(el);
                 Owners owners(self, v);
                 v.push_back(std::move(handle));
                 owners.append(el);
             },
             py::arg("element"),
             "Add a single element to the lattice.\n\n"
             "The element is shared, not copied: later changes to it are visible here, and\n"
             "appending the same element twice creates two occurrences of one element."
        )

        .def("extend",
             [](py::object self, py::iterable const & l) {
                 auto & v = self.cast<KnownElementsList &>();

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

        .def("size", &KnownElementsList::size)
        .def("is_empty", &KnownElementsList::empty)
        .def("__len__", [](const KnownElementsList &v) { return v.size(); },
             "The length of the list.")

        .def("clear",
             [](py::object self) {
                 auto & v = self.cast<KnownElementsList &>();
                 Owners owners(self, v);
                 v.clear();
                 owners.clear();
             },
             "Remove all elements from the lattice."
        )

        .def("pop_back",
             [](py::object self) {
                 auto & v = self.cast<KnownElementsList &>();
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
                 auto & v = self.cast<KnownElementsList &>();
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
                 auto & v = self.cast<KnownElementsList &>();
                 auto const i = checked_index(v, index);
                 Owners owners(self, v);
                 return owners.get(i);
             },
             py::arg("index"),
             "Return the element at a position. This is the element itself, not a copy."
        )

        .def("__setitem__",
             [checked_index](py::object self, py::ssize_t index, py::object el) {
                 auto & v = self.cast<KnownElementsList &>();
                 auto const i = checked_index(v, index);
                 auto handle = handle_from_python(el);
                 Owners owners(self, v);
                 v[i] = std::move(handle);
                 owners.set(i, el);
             },
             py::arg("index"), py::arg("element"),
             "Replace the element at a position."
        )

        .def("__getitem__",
             [](py::object self, py::slice const & slice) {
                 auto & v = self.cast<KnownElementsList &>();
                 size_t start = 0, stop = 0, step = 0, length = 0;
                 if (!slice.compute(v.size(), &start, &stop, &step, &length)) {
                     throw py::error_already_set();
                 }
                 Owners owners(self, v);

                 // a slice is a new lattice over the same elements, like a Python list slice
                 py::object out = py::type::of(self)();
                 py::list picked;
                 for (size_t i = 0; i < length; ++i) { picked.append(owners.get(start + i * step)); }
                 out.attr("extend")(picked);
                 return out;
             },
             py::arg("slice"),
             "Return a new lattice over the selected elements. The elements are shared."
        )

        .def("__setitem__",
             [](py::object self, py::slice const & slice, py::iterable const & value) {
                 auto & v = self.cast<KnownElementsList &>();
                 size_t start = 0, stop = 0, step = 0, length = 0;
                 if (!slice.compute(v.size(), &start, &stop, &step, &length)) {
                     throw py::error_already_set();
                 }

                 // convert first: a bad entry must leave the lattice untouched
                 std::vector<std::pair<elements::ElementHandle, py::object>> pending;
                 for (auto const & item : value)
                 {
                     py::object el = py::reinterpret_borrow<py::object>(item);
                     pending.emplace_back(handle_from_python(el), el);
                 }

                 if (step != 1 && pending.size() != length)
                 {
                     throw py::value_error(
                         "attempt to assign sequence of size " + std::to_string(pending.size()) +
                         " to extended slice of size " + std::to_string(length));
                 }

                 Owners owners(self, v);
                 if (step == 1)
                 {
                     // a contiguous slice may change the length
                     for (size_t i = 0; i < length; ++i) { v.erase(start); owners.erase(start); }
                     for (size_t i = 0; i < pending.size(); ++i)
                     {
                         v.insert(start + i, std::move(pending[i].first));
                         owners.insert(start + i, pending[i].second);
                     }
                 }
                 else
                 {
                     for (size_t i = 0; i < length; ++i)
                     {
                         auto const at = start + i * step;
                         v[at] = std::move(pending[i].first);
                         owners.set(at, pending[i].second);
                     }
                 }
             },
             py::arg("slice"), py::arg("elements"),
             "Replace the selected elements."
        )

        .def("__delitem__",
             [](py::object self, py::slice const & slice) {
                 auto & v = self.cast<KnownElementsList &>();
                 size_t start = 0, stop = 0, step = 0, length = 0;
                 if (!slice.compute(v.size(), &start, &stop, &step, &length)) {
                     throw py::error_already_set();
                 }
                 Owners owners(self, v);

                 // back to front, so earlier positions stay valid as we remove
                 for (size_t i = length; i-- > 0; )
                 {
                     auto const at = start + i * step;
                     v.erase(at);
                     owners.erase(at);
                 }
             },
             py::arg("slice"),
             "Remove the selected elements."
        )

        .def("__delitem__",
             [checked_index](py::object self, py::ssize_t index) {
                 auto & v = self.cast<KnownElementsList &>();
                 auto const i = checked_index(v, index);
                 Owners owners(self, v);
                 v.erase(i);
                 owners.erase(i);
             },
             py::arg("index"),
             "Remove the element at a position.\n\n"
             "If the same element occupies other positions, those are untouched: this\n"
             "removes one occurrence, not the element."
        )

        .def("insert",
             [](py::object self, py::ssize_t index, py::object el) {
                 auto & v = self.cast<KnownElementsList &>();
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

        // Identity first, then equality. Elements compare by value, so two distinct
        // elements with the same parameters are equal; for locating an occurrence in a
        // lattice, the object the caller means is the one they hold.
        .def("index",
             [](py::object self, py::object el) {
                 auto & v = self.cast<KnownElementsList &>();
                 Owners owners(self, v);
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                 {
                     if (owners.get(i).is(el)) { return i; }
                 }
                 throw py::value_error("element is not in the lattice");
             },
             py::arg("element"),
             "Return the first position holding this element."
        )

        .def("count",
             [](py::object self, py::object el) {
                 auto & v = self.cast<KnownElementsList &>();
                 Owners owners(self, v);
                 KnownElementsList::size_type n = 0;
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                 {
                     if (owners.get(i).is(el)) { ++n; }
                 }
                 return n;
             },
             py::arg("element"),
             "Return how many positions this element occupies."
        )

        .def("__contains__",
             [](py::object self, py::object el) {
                 auto & v = self.cast<KnownElementsList &>();
                 Owners owners(self, v);
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                 {
                     if (owners.get(i).is(el)) { return true; }
                 }
                 return false;
             },
             py::arg("element")
        )

        .def("remove",
             [](py::object self, py::object el) {
                 auto & v = self.cast<KnownElementsList &>();
                 Owners owners(self, v);
                 for (KnownElementsList::size_type i = 0; i < v.size(); ++i)
                 {
                     if (owners.get(i).is(el))
                     {
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
                 auto & v = self.cast<KnownElementsList &>();
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
                KnownElementsList const & v,
                RefPart ref,
                std::string order,
                bool fallback_identity_map
            )
            {
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
            [](KnownElementsList const & v, RefPart ref)  // intentional copy of ref
            {
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
        &impactx::elements::transformation::insert_element_every_ds,
        py::arg("list"),
        py::arg("ds"),
        py::arg("element"),
        "Insert an element every s into an element list"
    );
}
