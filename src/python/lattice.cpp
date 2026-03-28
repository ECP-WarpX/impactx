/* Copyright 2022-2026 The ImpactX Community
 *
 * Authors: Axel Huebl, Eric G. Stern
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <diagnostics/LinearMap.H>
#include <elements/All.H>
#include <elements/mixin/lineartransport.H>
#include <elements/transformation/Insert.H>
#include <particles/CovarianceMatrix.H>
#include <particles/ReferenceParticle.H>

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

    // all-element type list
    //
    // Ref semantics (Python vs this binding): a Python list keeps references to mutable
    // objects, but here each append/extend materializes a C++ KnownElements value in this
    // list (copy or move). The original Python handle is a separate object unless we
    // intentionally bind shared ownership (e.g. std::shared_ptr element types later).
    //
    // TODO: register GPU-heavy element types with py::class_<T, std::shared_ptr<T>> and
    //       store std::shared_ptr in KnownElements (or a thin wrapper) so append/extend
    //       match Python mutable-object reference semantics.
    // TODO: optional weak_ptr to a SimulationLifetime/session object on elements for
    //       clear errors on use after sim.finalize() (coordinate with pyAMReX session work).
    //       The shared token should be created with ImpactX construction (not init_grids),
    //       and expose readiness (amrex_ready) so pre-init objects can exist but reject
    //       GPU use until sim.init_grids() makes AMReX ready.
    using KnownElementsList = std::list<KnownElements>;
    py::class_<KnownElementsList> kel(me, "KnownElementsList");
    kel
        .def(py::init<>())
        .def(py::init<KnownElements>())
        .def(py::init([](py::list const & l){
            KnownElementsList v;
            for (auto const & handle : l)
                v.push_back(handle.cast<KnownElements>());
            return v;  // return by value
        }))

        .def("append", [](KnownElementsList &v, KnownElements el) { v.emplace_back(std::move(el)); },
             "Add a single element to the list.")

        .def("extend",
             [](KnownElementsList &v, KnownElementsList const & l) {
                 for (auto const & el : l)
                     v.push_back(el);
                 return v;
             },
             "Add a list of elements to the list.")
        .def("extend",
             [](KnownElementsList &v, py::list const & l) {
                 for (auto const & handle : l)
                 {
                     auto el = handle.cast<KnownElements>();
                     v.push_back(el);
                 }
                 return v;
             },
             "Add a list of elements to the list."
        )

        .def("size", &KnownElementsList::size)
        .def("clear", &KnownElementsList::clear,
             "Clear the list to become empty.")
        .def("is_empty", &KnownElementsList::empty)
        .def("pop_back", &KnownElementsList::pop_back,
             "Return and remove the last element of the list.")
        .def("__len__", [](const KnownElementsList &v) { return v.size(); },
             "The length of the list.")
        .def("__iter__", [](KnownElementsList &v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())  // Keep list alive while iterator is used
        .def("__getitem__", [](KnownElementsList &v, size_t index) -> elements::KnownElements& {
            if (index >= v.size()) {
                throw std::out_of_range("Index out of range");
            }
            auto it = std::next(v.begin(), index);
            return *it;  // return by reference
        }, py::return_value_policy::reference_internal)

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
