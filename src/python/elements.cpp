/* Copyright 2022-2026 The ImpactX Community
 *
 * Authors: Axel Huebl, Eric G. Stern
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/Push.H>
#include <elements/All.H>
#include <elements/mixin/accessors.H>
#include <particles/CovarianceMatrix.H>

#include <AMReX_Enum.H>
#include <AMReX_REAL.H>

#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace py = pybind11;
using namespace impactx;


namespace detail
{
    /** Helper Function for Property Getters
     *
     * This queries an amrex::ParmParse entry. This throws a
     * std::runtime_error if the entry is not found.
     *
     * This handles most common throw exception logic in ImpactX instead of
     * going over library boundaries via amrex::Abort().
     *
     * @tparam T type of the amrex::ParmParse entry
     * @param prefix the prefix, e.g., "impactx" or "amr"
     * @param name the actual key of the entry, e.g., "particle_shape"
     * @return the queried value (or throws if not found)
     */
    template< typename T>
    auto get_or_throw (std::string const & prefix, std::string const & name)
    {
        using V = std::decay_t<T>;
        V value;

        bool has_name = false;
        // TODO: if array do queryarr
        // has_name = amrex::ParmParse(prefix).queryarr(name.c_str(), value);
        if constexpr (std::is_same_v<V, bool> || std::is_same_v<V, std::string>)
            has_name = amrex::ParmParse(prefix).query(name.c_str(), value);
        else
            has_name = amrex::ParmParse(prefix).queryWithParser(name.c_str(), value);

        if (!has_name)
            throw std::runtime_error(prefix + "." + name + " is not set yet");
        return value;
    }
}

namespace
{
    /** Registers the mixin BeamOptics::operator methods
     */
    template<typename T_PyClass>
    void register_beamoptics_push (T_PyClass & cl)
    {
        using Element = typename T_PyClass::type;  // py::class<T, options...>

        cl.def("push",
            [](Element & el, ImpactXParticleContainer & pc, int step, int period) {
                el(pc, step, period);
            },
            py::arg("pc"), py::arg("step")=0, py::arg("period")=0,
            "Push first the reference particle, then all other particles."
        )
        .def("finalize", &Element::finalize)
        ;
    }

    /** Registers the mixin LinearTransport::operator method
     */
    template<typename T_PyClass>
    void register_envelope_push (T_PyClass & cl)
    {
        using Element = typename T_PyClass::type;  // py::class<T, options...>

        cl.def("push",
               [](Element & el, Map6x6 & cm, RefPart const & ref) {
                   el(cm, ref);
               },
               py::arg("cm"), py::arg("ref"),
               "Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first."
        );
    }

    /** Register the transfer_map() method
     *
     * Exposes the element's own analytic linear transport map (see the
     * @c mixin::LinearTransport CRTP mixin). An element that implements a linear map
     * (@c has_linear_transport) returns it, one that does not throws a uniform,
     * self-documenting error.
     *
     * This is the element's per-slice map (for @c ds/nslice), i.e. the building
     * block that @c KnownElementsList.transfer_map composes over @c nslice and
     * over the whole lattice.
     *
     * The evaluation convention matches @c diagnostics::walk_lattice and the
     * envelope tracker: advance a private copy of the reference particle by one
     * slice, then evaluate @c transport_map at the post-advance reference state.
     * Elements whose linear map depends on the reference momentum through the
     * slice (@c ShortRF, @c ChrAcc) read the post-advance @c pt, and elements
     * that integrate their map during the reference push (@c RFCavity,
     * @c SoftQuadrupole, @c SoftSolenoid) only have one to return afterwards.
     */
    template<typename T_PyClass>
    void register_transfer_map (T_PyClass & cl)
    {
        using Element = typename T_PyClass::type;  // py::class<T, options...>

        cl.def("transfer_map",
            [](Element const & el, RefPart const & ref) -> Map6x6 {
                if constexpr (Element::has_linear_transport) {
                    // private copy: the caller's reference particle is not modified
                    RefPart ref_slice = ref;
                    // the element starts at the current s, as in walk_lattice
                    ref_slice.set_edge();
                    el(ref_slice);
                    return el.transport_map(ref_slice);
                } else {
                    throw std::runtime_error(
                        std::string(Element::type)
                        + ": Linear transport map is not yet implemented for this element."
                    );
                }
            },
            py::arg("ref"),
            "Return this element's 6x6 linear transport map for the given\n"
            "reference particle.\n\n"
            "Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).\n"
            "For an element with ``nslice`` > 1 this is the map of a single\n"
            "``ds/nslice`` slice (the building block that the lattice transfer\n"
            "map composes). Raises for an element whose linear transport map is\n"
            "not implemented.\n\n"
            "Uses the same convention as ``KnownElementsList.transfer_map``: the\n"
            "reference particle is advanced through one slice first, and the map\n"
            "is evaluated at that post-advance reference state. The reference\n"
            "particle you pass in is copied and left unmodified.\n\n"
            ":param ref: reference particle at the element entrance\n"
        );
    }

    /** Register push() method overloads */
    template<typename T_PyClass>
    void register_push (T_PyClass & cl)
    {
        register_beamoptics_push(cl);
        register_envelope_push(cl);
        register_transfer_map(cl);
    }

    /** Register reverse() method */
    template<typename T_PyClass>
    void register_reverse (T_PyClass & cl)
    {
        using Element = typename T_PyClass::type;

        cl.def("reverse", &Element::reverse,
            "Reverse the element in-place so that pushing particles through\n"
            "it reverses the effect of the original element."
        );
    }

    /** Trait for std::optional values, which are formatted as "None" if unset */
    template<typename T>
    struct is_optional : std::false_type {};

    template<typename T>
    struct is_optional<std::optional<T>> : std::true_type {};

    /** Helper to format {key, value} pairs
     *
     * Expected outcome is ", key=value" with key as a string and appropriate formatting for value;
     * an unset std::optional value is formatted as "None".
     *
     * @tparam T value type
     * @param arg a key-value pair
     * @return a string of the form ", key=value"
     */
    template<typename T>
    std::string
    format_extra (std::pair<char const *, T> const & arg)
    {
        if constexpr (is_optional<T>::value)
        {
            // optional: unset is None, otherwise the value it holds
            if (!arg.second.has_value()) {
                return std::string(", ")
                    .append(arg.first)
                    .append("=None");
            }
            return format_extra(std::make_pair(arg.first, arg.second.value()));
        }
        else if constexpr (std::is_floating_point_v<T>)
        {
            // float
            // TODO: format as scientific number
            return std::string(", ")
                .append(arg.first)
                .append("=")
                .append(std::to_string(arg.second));
        }
        else if constexpr (std::is_integral_v<T>)
        {
            // int
            return std::string(", ")
                .append(arg.first)
                .append("=")
                .append(std::to_string(arg.second));
        } else
        {
            // already a string
            return std::string(", ")
                .append(arg.first)
                .append("=")
                .append(arg.second);
        }
    }

    /** Helper to build a __repr__ for an Element
     *
     * @tparam T_Element type for the C++ element type
     * @tparam ExtraArgs type for pairs of name, value to add
     * @param el the current element
     * @param args pars of name, value to add
     * @return a string suitable for Python's __repr__
     */
    template<typename T_Element, typename... ExtraArgs>
    std::string element_name (T_Element const & el, std::pair<char const *, ExtraArgs> const &... args)
    {
        // Fixed element type name, e.g., "SBend"
        std::string const type = T_Element::type;

        // User-provided element name, e.g., "name=bend1"
        std::string const name = el.has_name() ? ", name=" + el.name() : "";

        // Noteworthy element parameters, e.g., "ds=2.3, key=value, ..."
        std::string extra_args;

        // properties of mixin classes
        if constexpr (elements::mixin::is_thick_v<T_Element>) {
            extra_args.append(format_extra(std::make_pair("ds", el.ds())));
        }

        // select properties specific to the element
        ((extra_args.append(format_extra(args))), ...);

        // combine it all together
        return "<impactx.elements." +
               type +
               name +
               extra_args +
               ">";
    }

    // List of property types we store in elements
    using ElementPropertyTypes = std::variant<
        amrex::ParticleReal,
        int,
        long,
        bool,
        std::string,
        std::vector<amrex::ParticleReal>,
        std::vector<int>,
        std::vector<long>,
        Map6x6,
        Vector3,
        Map3x6,
        py::none
    >;

    /** Create a map (for a Python dictionary) of properties of an element
     *
     * The values are meant to align with the element type's Python constructor keyword
     * arguments where possible, used internally for copies and (de)serialization.
     *
     * The dict always includes a ``type`` string (the element class name) for dispatch;
     * we also include ``ds`` as 0.0 for thin elements for simplicity (plots, etc.);
     * ``type`` is not a constructor argument and neither is ``ds`` for thin elements,
     * and must be omitted when unpacking, e.g.:
     * ```py
     * dr = elements.Drift(name="drift1", ds=1.0)
     * d = dr.to_dict()
     * kwargs = {k: v for k, v in d.items() if k != "type" and (k != "ds" or v != 0.0)}
     * dr2 = elements.Drift(**kwargs)
     * ```
     */
    template<typename T_Element, typename... ExtraArgs>
    std::map<std::string, ElementPropertyTypes>
    element_dict (T_Element const & el, std::pair<char const *, ExtraArgs> const &... args)
    {
        // Fixed element type name, e.g., "SBend"
        std::string const type = T_Element::type;

        ElementPropertyTypes name = py::none();
        if (el.has_name()) { name = el.name(); }

        // combine it all together
        std::map<std::string, ElementPropertyTypes> ret = {
            {"type", type},
            {"name", name},
            // Thin / Thick properties
            {"ds", el.ds()}
        };

        // properties of mixin classes
        if (el.nslice() != 1) {
            // for thick elements, nslice default is 1 and can be omitted
            // for thin elements, nslice is 1 by definition and not a constructor argument
            ret.insert(std::make_pair("nslice", el.nslice()));
        }
        if constexpr (elements::mixin::has_alignment_v<T_Element>) {
            ret.insert(std::make_pair("dx", el.dx()));
            ret.insert(std::make_pair("dy", el.dy()));
            ret.insert(std::make_pair("rotation", el.rotation()));
        }
        if constexpr (elements::mixin::has_pipe_aperture_v<T_Element>) {
            ret.insert(std::make_pair("aperture_x", el.aperture_x()));
            ret.insert(std::make_pair("aperture_y", el.aperture_y()));
        }

        // properties specific to the element
        ((ret.insert(args)), ...);

        return ret;
    }
}


void init_lattice(py::module&);

void init_elements(py::module& m)
{
    using namespace elements;

    m.attr("Map6x6") = py::type::of<Map6x6>();
    m.attr("Map3x6") = py::type::of<Map3x6>();
    m.attr("Vector3") = py::type::of<Vector3>();

    py::module_ me = m.def_submodule(
        "elements",
        "Accelerator lattice elements in ImpactX"
    );

    // mixin classes

    py::module_ const mx = me.def_submodule(
        "mixin",
        "Mixin classes for accelerator lattice elements in ImpactX"
    );

    py::class_<elements::mixin::Named>(mx, "Named")
        .def_property("name",
            [](elements::mixin::Named & nm) -> std::optional<std::string> {
                return nm.has_name() ? std::optional<std::string>{nm.name()} : std::nullopt;
            },
            [](elements::mixin::Named & nm, std::string new_name) { nm.set_name(new_name); },
            "segment length in m"
        )
        .def_property_readonly("has_name", &elements::mixin::Named::has_name)
    ;

    py::class_<elements::mixin::Thick>(mx, "Thick")
        .def_property("ds",
            [](elements::mixin::Thick & th) { return th.m_ds; },
            [](elements::mixin::Thick & th, amrex::ParticleReal ds) { th.m_ds = ds; },
            "segment length in m"
        )
        .def_property("nslice",
            [](elements::mixin::Thick & th) { return th.m_nslice; },
            [](elements::mixin::Thick & th, int nslice) {
                if (nslice <= 0)
                    throw std::invalid_argument("Thick: nslice must be > 0!");
                th.m_nslice = nslice;
            },
            "number of slices used for the application of space charge"
        )
    ;

    py::class_<elements::mixin::Thin>(mx, "Thin")
        .def_property_readonly("ds",
            &elements::mixin::Thin::ds,
            "segment length in m"
        )
        .def_property_readonly("nslice",
            &elements::mixin::Thin::nslice,
            "number of slices used for the application of space charge"
        )
    ;

    py::class_<elements::mixin::Alignment>(mx, "Alignment")
        .def_property("dx",
            [](elements::mixin::Alignment & a) { return a.dx(); },
            [](elements::mixin::Alignment & a, amrex::ParticleReal dx) { a.m_dx = dx; },
            "horizontal translation error in m"
        )
        .def_property("dy",
            [](elements::mixin::Alignment & a) { return a.dy(); },
            [](elements::mixin::Alignment & a, amrex::ParticleReal dy) { a.m_dy = dy; },
            "vertical translation error in m"
        )
        .def_property("rotation",
            [](elements::mixin::Alignment & a) { return a.rotation(); },
            [](elements::mixin::Alignment & a, amrex::ParticleReal rotation_degree)
            {
                a.m_rotation = rotation_degree * elements::mixin::Alignment::degree2rad;
            },
            "rotation error in the transverse plane in degree"
        )
    ;

    py::class_<elements::mixin::PipeAperture>(mx, "PipeAperture")
        .def_property("aperture_x",
            [](elements::mixin::PipeAperture & pa) { return pa.aperture_x(); },
            [](elements::mixin::PipeAperture & pa, amrex::ParticleReal aperture_x)
            {
                pa.set_aperture_x(aperture_x);
            },
            "horizontal aperture in m; zero or less removes the aperture restriction"
        )
        .def_property("aperture_y",
            [](elements::mixin::PipeAperture & pa) { return pa.aperture_y(); },
            [](elements::mixin::PipeAperture & pa, amrex::ParticleReal aperture_y)
            {
                pa.set_aperture_y(aperture_y);
            },
            "vertical aperture in m; zero or less removes the aperture restriction"
        )
    ;

    /* TODO
    py::class_<elements::mixin::LinearTransport>(mx, "LinearTransport")
        // type of map
        .def_property_readonly_static("Map6x6",
              [](py::object){ return py::type::of<elements::mixin::LinearTransport::R>(); },
              "1-indexed, Fortran-ordered, 6x6 linear transport map type"
        )
        // values of the map
        //.def_property_readonly("R",
        //      [](elements::mixin::LinearTransport const & lt) { return lt.m_transport_map; },
        //      "1-indexed, Fortran-ordered, 6x6 linear transport map values"
        //)
    ;
    */

    // diagnostics

    py::class_<diagnostics::BeamMonitor, elements::mixin::Thin> py_BeamMonitor(me, "BeamMonitor");
    py_BeamMonitor
        .def("__repr__",
             [](diagnostics::BeamMonitor const & bm) {
                 return element_name(bm);
             }
        )
        .def("to_dict",
             [](diagnostics::BeamMonitor const & bm) {
                 return element_dict(
                     bm,
                     std::make_pair("backend", bm.backend()),
                     std::make_pair("encoding", bm.encoding()),
                     std::make_pair("period_sample_intervals", bm.period_sample_intervals()),
                     std::make_pair("particles", bm.particles()),
                     std::make_pair("beam_moments", bm.beam_moments())
                 );
             }
        )
        .def(py::init<std::string, std::string, std::string, int, bool, bool>(),
             py::arg("name"),
             py::arg("backend") = diagnostics::BeamMonitor::DEFAULT_backend,
             py::arg("encoding") = diagnostics::BeamMonitor::DEFAULT_encoding,
             py::arg("period_sample_intervals") =
                 diagnostics::BeamMonitor::DEFAULT_period_sample_intervals,
             py::arg("particles") = diagnostics::BeamMonitor::DEFAULT_particles,
             py::arg("beam_moments") = diagnostics::BeamMonitor::DEFAULT_beam_moments,
             "This element writes the particle beam out to openPMD data."
        )
        .def_property_readonly("name",
            &diagnostics::BeamMonitor::name,
            "name of the series"
        )
        .def_property_readonly("has_name", &diagnostics::BeamMonitor::has_name)
        .def_property_readonly("backend",
            &diagnostics::BeamMonitor::backend,
            "openPMD file backend (e.g. default, bp4, h5)"
        )
        .def_property_readonly("encoding",
            &diagnostics::BeamMonitor::encoding,
            R"(openPMD iteration encoding: "v" variable-based, "f" file-based, "g" group-based)"
        )
        .def_property_readonly("period_sample_intervals",
            &diagnostics::BeamMonitor::period_sample_intervals,
            "for periodic lattices, only output every Nth period (turn or cycle)"
        )
        .def_property("particles",
            [](diagnostics::BeamMonitor & bm) { return bm.particles(); },
            [](diagnostics::BeamMonitor & bm, bool particles) { bm.set_particles(particles); },
            "Output the individual particles of the beam (default: True)"
        )
        .def_property("beam_moments",
            [](diagnostics::BeamMonitor & bm) { return bm.beam_moments(); },
            [](diagnostics::BeamMonitor & bm, bool beam_moments) { bm.set_beam_moments(beam_moments); },
            "Output the beam moments (reduced beam characteristics) as openPMD attributes (default: True)"
        )
        .def_property("nonlinear_lens_invariants",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<bool>(bm.name(), "nonlinear_lens_invariants"); },
            [](diagnostics::BeamMonitor & bm, bool nonlinear_lens_invariants) {
                amrex::ParmParse pp_element(bm.name());
                pp_element.add("nonlinear_lens_invariants", nonlinear_lens_invariants);
            },
            "Compute and output the invariants H and I within the nonlinear magnetic insert element. "
            "Requires particles=True."
        )
        .def_property("alpha",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.name(), "alpha"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real alpha) {
                amrex::ParmParse pp_element(bm.name());
                pp_element.add("alpha", alpha);
            },
            "Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.\n"
            "Horizontal and vertical values must be equal."
        )
        .def_property("beta",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.name(), "beta"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real beta) {
                amrex::ParmParse pp_element(bm.name());
                pp_element.add("beta", beta);
            },
            "Twiss beta (in meters) of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.\n"
            "Horizontal and vertical values must be equal."
        )
        .def_property("tn",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.name(), "tn"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real tn) {
                amrex::ParmParse pp_element(bm.name());
                pp_element.add("tn", tn);
            },
            "Dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I."
        )
        .def_property("cn",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.name(), "cn"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real cn) {
                amrex::ParmParse pp_element(bm.name());
                pp_element.add("cn", cn);
            },
            "Scale factor (in meters^(1/2)) of the IOTA nonlinear magnetic insert element used for computing H and I."
        )
    ;
    register_push(py_BeamMonitor);
    register_reverse(py_BeamMonitor);

    // beam optics

    py::class_<Aperture, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment> py_Aperture(me, "Aperture");
    py_Aperture
        .def("__repr__",
             [](Aperture const & ap) {
                 return element_name(
                    ap,
                    std::make_pair("shape", amrex::getEnumNameString(ap.m_shape)),
                    std::make_pair("action", amrex::getEnumNameString(ap.m_action))
                );
             }
        )
        .def("to_dict",
            [](Aperture const & ap) {
                return element_dict(
                    ap,
                    std::make_pair("shape", amrex::getEnumNameString(ap.m_shape)),
                    std::make_pair("action", amrex::getEnumNameString(ap.m_action)),
                    std::make_pair("aperture_x", ap.aperture_x()),
                    std::make_pair("aperture_y", ap.aperture_y()),
                    std::make_pair("repeat_x", ap.m_repeat_x),
                    std::make_pair("repeat_y", ap.m_repeat_y),
                    std::make_pair("shift_odd_x", ap.m_shift_odd_x)
                );
            }
        )
        .def(py::init([](
                 amrex::ParticleReal aperture_x,
                 amrex::ParticleReal aperture_y,
                 amrex::ParticleReal repeat_x,
                 amrex::ParticleReal repeat_y,
                 bool shift_odd_x,
                 std::string const & shape,
                 std::string const & action,
                 amrex::ParticleReal dx,
                 amrex::ParticleReal dy,
                 amrex::ParticleReal rotation_degree,
                 std::optional<std::string> name
             )
             {
                 Aperture::Shape const s = amrex::getEnum<Aperture::Shape>(shape);
                 Aperture::Action const a = amrex::getEnum<Aperture::Action>(action);
                 return new Aperture(aperture_x, aperture_y, repeat_x, repeat_y, shift_odd_x, s, a, dx, dy, rotation_degree, name);
             }),
             py::arg("aperture_x"),
             py::arg("aperture_y"),
             py::arg("repeat_x") = Aperture::DEFAULT_repeat_x,
             py::arg("repeat_y") = Aperture::DEFAULT_repeat_y,
             py::arg("shift_odd_x") = Aperture::DEFAULT_shift_odd_x,
             py::arg("shape") = amrex::getEnumNameString(Aperture::DEFAULT_shape),
             py::arg("action") = amrex::getEnumNameString(Aperture::DEFAULT_action),
             py::arg("dx") = Aperture::DEFAULT_dx,
             py::arg("dy") = Aperture::DEFAULT_dy,
             py::arg("rotation") = Aperture::DEFAULT_rotation_degree,
             py::arg("name") = py::none(),
             "A short collimator element applying a transverse aperture boundary."
        )
        .def_property("shape",
            [](Aperture & ap)
            {
                return amrex::getEnumNameString(ap.m_shape);
            },
            [](Aperture & ap, std::string const & shape)
            {
                ap.m_shape = amrex::getEnum<Aperture::Shape>(shape);
            },
            "aperture type (rectangular, elliptical)"
        )
        .def_property("action",
            [](Aperture & ap)
            {
                return amrex::getEnumNameString(ap.m_action);
            },
            [](Aperture & ap, std::string const & action)
            {
                ap.m_action = amrex::getEnum<Aperture::Action>(action);
            },
            "action type (transmit, absorb)"
        )
        .def_property("aperture_x",
            [](Aperture & ap) { return ap.aperture_x(); },
            [](Aperture & ap, amrex::ParticleReal aperture_x) { ap.set_aperture_x(aperture_x); },
            "maximum horizontal coordinate in m; zero or less removes the horizontal boundary"
        )
        .def_property("aperture_y",
            [](Aperture & ap) { return ap.aperture_y(); },
            [](Aperture & ap, amrex::ParticleReal aperture_y) { ap.set_aperture_y(aperture_y); },
            "maximum vertical coordinate in m; zero or less removes the vertical boundary"
        )
        .def_property("repeat_x",
            [](Aperture & ap) { return ap.m_repeat_x; },
            [](Aperture & ap, amrex::ParticleReal repeat_x) { ap.m_repeat_x = repeat_x; },
            "horizontal period for repeated aperture masking"
        )
        .def_property("repeat_y",
            [](Aperture & ap) { return ap.m_repeat_y; },
            [](Aperture & ap, amrex::ParticleReal repeat_y) { ap.m_repeat_y = repeat_y; },
            "vertical period for repeated aperture masking"
        )
        .def_property("shift_odd_x",
            [](Aperture & ap) { return ap.m_shift_odd_x; },
            [](Aperture & ap, bool shift_odd_x) { ap.m_shift_odd_x = shift_odd_x; },
            "for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. "
            "Use alignment offsets dx,dy to move whole mask as needed."
        )
    ;
    register_push(py_Aperture);
    register_reverse(py_Aperture);

    py::class_<ChrDrift, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ChrDrift(me, "ChrDrift");
    py_ChrDrift
        .def("__repr__",
             [](ChrDrift const & chr_drift) {
                 return element_name(chr_drift);
             }
        )
        .def("to_dict",
            [](ChrDrift const & chr_drift) {
                return element_dict(chr_drift);
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("dx") = ChrDrift::DEFAULT_dx,
             py::arg("dy") = ChrDrift::DEFAULT_dy,
             py::arg("rotation") = ChrDrift::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ChrDrift::DEFAULT_aperture_x,
             py::arg("aperture_y") = ChrDrift::DEFAULT_aperture_y,
             py::arg("nslice") = ChrDrift::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A Drift with chromatic effects included."
        )
    ;
    register_push(py_ChrDrift);
    register_reverse(py_ChrDrift);

    py::class_<ChrQuad, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ChrQuad(me, "ChrQuad");
    py_ChrQuad
        .def("__repr__",
             [](ChrQuad const & chr_quad) {
                 return element_name(
                     chr_quad,
                     std::make_pair("k", chr_quad.m_k)
                 );
             }
        )
        .def("to_dict",
             [](ChrQuad const & chr_quad) {
                 return element_dict(
                     chr_quad,
                     std::make_pair("k", chr_quad.m_k),
                     std::make_pair("unit", chr_quad.m_unit)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("k"),
             py::arg("unit") = ChrQuad::DEFAULT_unit,
             py::arg("dx") = ChrQuad::DEFAULT_dx,
             py::arg("dy") = ChrQuad::DEFAULT_dy,
             py::arg("rotation") = ChrQuad::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ChrQuad::DEFAULT_aperture_x,
             py::arg("aperture_y") = ChrQuad::DEFAULT_aperture_y,
             py::arg("nslice") = ChrQuad::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A Quadrupole magnet with chromatic effects included."
        )
        .def_property("k",
            [](ChrQuad & cq) { return cq.m_k; },
            [](ChrQuad & cq, amrex::ParticleReal k) { cq.m_k = k; },
            "quadrupole strength in 1/m^2 (or T/m)"
        )
        .def_property("unit",
            [](ChrQuad & cq) { return cq.m_unit; },
            [](ChrQuad & cq, int unit) { cq.m_unit = unit; },
            "unit specification for quad strength"
        )
    ;
    register_push(py_ChrQuad);
    register_reverse(py_ChrQuad);

    py::class_<ChrPlasmaLens, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ChrPlasmaLens(me, "ChrPlasmaLens");
    py_ChrPlasmaLens
        .def("__repr__",
             [](ChrPlasmaLens const & chr_pl_lens) {
                 return element_name(
                     chr_pl_lens,
                     std::make_pair("k", chr_pl_lens.m_k)
                 );
             }
        )
        .def("to_dict",
             [](ChrPlasmaLens const & chr_pl_lens) {
                 return element_dict(
                     chr_pl_lens,
                     std::make_pair("k", chr_pl_lens.m_k),
                     std::make_pair("unit", chr_pl_lens.m_unit)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("k"),
             py::arg("unit") = ChrPlasmaLens::DEFAULT_unit,
             py::arg("dx") = ChrPlasmaLens::DEFAULT_dx,
             py::arg("dy") = ChrPlasmaLens::DEFAULT_dy,
             py::arg("rotation") = ChrPlasmaLens::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ChrPlasmaLens::DEFAULT_aperture_x,
             py::arg("aperture_y") = ChrPlasmaLens::DEFAULT_aperture_y,
             py::arg("nslice") = ChrPlasmaLens::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "An active Plasma Lens with chromatic effects included."
        )
        .def_property("k",
            [](ChrPlasmaLens & chr_pl_lens) { return chr_pl_lens.m_k; },
            [](ChrPlasmaLens & chr_pl_lens, amrex::ParticleReal k) { chr_pl_lens.m_k = k; },
            "focusing strength in 1/m^2 (or T/m)"
        )
        .def_property("unit",
            [](ChrPlasmaLens & chr_pl_lens) { return chr_pl_lens.m_unit; },
            [](ChrPlasmaLens & chr_pl_lens, int unit) { chr_pl_lens.m_unit = unit; },
            "unit specification for focusing strength"
        )
    ;
    register_push(py_ChrPlasmaLens);
    register_reverse(py_ChrPlasmaLens);

    py::class_<ChrAcc, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment> py_ChrAcc(me, "ChrAcc");
    py_ChrAcc
        .def("__repr__",
             [](ChrAcc const & chr_acc) {
                 return element_name(
                     chr_acc,
                     std::make_pair("ez", chr_acc.m_ez),
                     std::make_pair("bz", chr_acc.m_bz)
                 );
             }
        )
        .def("to_dict",
             [](ChrAcc const & chr_acc) {
                 return element_dict(
                     chr_acc,
                     std::make_pair("ez", chr_acc.m_ez),
                     std::make_pair("bz", chr_acc.m_bz)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("ez"),
             py::arg("bz"),
             py::arg("dx") = ChrAcc::DEFAULT_dx,
             py::arg("dy") = ChrAcc::DEFAULT_dy,
             py::arg("rotation") = ChrAcc::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ChrAcc::DEFAULT_aperture_x,
             py::arg("aperture_y") = ChrAcc::DEFAULT_aperture_y,
             py::arg("nslice") = ChrAcc::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A region of Uniform Acceleration, with chromatic effects included."
        )
        .def_property("ez",
            [](ChrAcc & ca) { return ca.m_ez; },
            [](ChrAcc & ca, amrex::ParticleReal ez) { ca.m_ez = ez; },
            "electric field strength in 1/m"
        )
        .def_property("bz",
            [](ChrAcc & ca) { return ca.m_bz; },
            [](ChrAcc & ca, amrex::ParticleReal bz) { ca.m_bz = bz; },
            "magnetic field strength in 1/m"
        )
    ;
    register_push(py_ChrAcc);
    register_reverse(py_ChrAcc);

    py::class_<ConstF, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ConstF(me, "ConstF");
    py_ConstF
        .def("__repr__",
             [](ConstF const & constf) {
                 return element_name(
                     constf,
                     std::make_pair("kx", constf.m_kx),
                     std::make_pair("ky", constf.m_ky),
                     std::make_pair("kt", constf.m_kt)
                 );
             }
        )
        .def("to_dict",
             [](ConstF const & constf) {
                 return element_dict(
                     constf,
                     std::make_pair("kx", constf.m_kx),
                     std::make_pair("ky", constf.m_ky),
                     std::make_pair("kt", constf.m_kt)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("kx"),
             py::arg("ky"),
             py::arg("kt"),
             py::arg("dx") = ConstF::DEFAULT_dx,
             py::arg("dy") = ConstF::DEFAULT_dy,
             py::arg("rotation") = ConstF::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ConstF::DEFAULT_aperture_x,
             py::arg("aperture_y") = ConstF::DEFAULT_aperture_y,
             py::arg("nslice") = ConstF::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A linear Constant Focusing element."
        )
        .def_property("kx",
            [](ConstF & cf) { return cf.m_kx; },
            [](ConstF & cf, amrex::ParticleReal kx) { cf.m_kx = kx; },
            "focusing x strength in 1/m"
        )
        .def_property("ky",
            [](ConstF & cf) { return cf.m_ky; },
            [](ConstF & cf, amrex::ParticleReal ky) { cf.m_ky = ky; },
            "focusing y strength in 1/m"
        )
        .def_property("kt",
            [](ConstF & cf) { return cf.m_kt; },
            [](ConstF & cf, amrex::ParticleReal kt) { cf.m_kt = kt; },
            "focusing t strength in 1/m"
        )
    ;
    register_push(py_ConstF);
    register_reverse(py_ConstF);

    py::class_<DipEdge, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_DipEdge(me, "DipEdge");
    py_DipEdge
        .def("__repr__",
             [](DipEdge const & dip_edge) {
                 return element_name(
                     dip_edge,
                     std::make_pair("psi", dip_edge.m_psi),
                     std::make_pair("rc", dip_edge.m_rc),
                     std::make_pair("g", dip_edge.m_g),
                     std::make_pair("R", dip_edge.m_R),
                     std::make_pair("K0", dip_edge.m_K0),
                     std::make_pair("K1", dip_edge.m_K1),
                     std::make_pair("K2", dip_edge.m_K2),
                     std::make_pair("K3", dip_edge.m_K3),
                     std::make_pair("K4", dip_edge.m_K4),
                     std::make_pair("K5", dip_edge.m_K5),
                     std::make_pair("K6", dip_edge.m_K6),
                     std::make_pair("model", amrex::getEnumNameString(dip_edge.m_model)),
                     std::make_pair("location", amrex::getEnumNameString(dip_edge.m_location)),
                     std::make_pair("modify_ref_part", dip_edge.m_modify_ref_part)
                 );
             }
        )
        .def("to_dict",
             [](DipEdge const & dip_edge) {
                 return element_dict(
                     dip_edge,
                     std::make_pair("psi", dip_edge.m_psi),
                     std::make_pair("rc", dip_edge.m_rc),
                     std::make_pair("g", dip_edge.m_g),
                     std::make_pair("R", dip_edge.m_R),
                     std::make_pair("K0", dip_edge.m_K0),
                     std::make_pair("K1", dip_edge.m_K1),
                     std::make_pair("K2", dip_edge.m_K2),
                     std::make_pair("K3", dip_edge.m_K3),
                     std::make_pair("K4", dip_edge.m_K4),
                     std::make_pair("K5", dip_edge.m_K5),
                     std::make_pair("K6", dip_edge.m_K6),
                     std::make_pair("model", amrex::getEnumNameString(dip_edge.m_model)),
                     std::make_pair("location", amrex::getEnumNameString(dip_edge.m_location)),
                     std::make_pair("modify_ref_part", dip_edge.m_modify_ref_part)
                 );
             }
        )
        .def(py::init([](
            amrex::ParticleReal psi,
            amrex::ParticleReal rc,
            amrex::ParticleReal g,
            amrex::ParticleReal R,
            amrex::ParticleReal K0,
            amrex::ParticleReal K1,
            amrex::ParticleReal K2,
            amrex::ParticleReal K3,
            amrex::ParticleReal K4,
            amrex::ParticleReal K5,
            amrex::ParticleReal K6,
            std::string const & model,
            std::string const & location,
            bool modify_ref_part,
            amrex::ParticleReal dx,
            amrex::ParticleReal dy,
            amrex::ParticleReal rotation_degree,
            amrex::ParticleReal aperture_x,
            amrex::ParticleReal aperture_y,
            std::optional<std::string> name
            )
            {
                 if (R <= 0.0)
                     throw std::invalid_argument(R"(DipEdge parameter R must be > 0.)");

                DipEdge::Model const fm = amrex::getEnum<DipEdge::Model>(model);
                DipEdge::Location const fl = amrex::getEnum<DipEdge::Location>(location);
                return new DipEdge(psi, rc, g, R, K0, K1, K2, K3, K4, K5, K6, fm, fl, modify_ref_part, dx, dy, rotation_degree, aperture_x, aperture_y, name);
            }),
            py::arg("psi"),
            py::arg("rc"),
            py::arg("g"),
            py::arg("R") = DipEdge::DEFAULT_R,
            py::arg("K0") = DipEdge::DEFAULT_K0,
            py::arg("K1") = DipEdge::DEFAULT_K1,
            py::arg("K2") = DipEdge::DEFAULT_K2,
            py::arg("K3") = DipEdge::DEFAULT_K3,
            py::arg("K4") = DipEdge::DEFAULT_K4,
            py::arg("K5") = DipEdge::DEFAULT_K5,
            py::arg("K6") = DipEdge::DEFAULT_K6,
            py::arg("model") = amrex::getEnumNameString(DipEdge::DEFAULT_model),
            py::arg("location") = amrex::getEnumNameString(DipEdge::DEFAULT_location),
            py::arg("modify_ref_part") = DipEdge::DEFAULT_modify_ref_part,
            py::arg("dx") = DipEdge::DEFAULT_dx,
            py::arg("dy") = DipEdge::DEFAULT_dy,
            py::arg("rotation") = DipEdge::DEFAULT_rotation_degree,
            py::arg("aperture_x") = DipEdge::DEFAULT_aperture_x,
            py::arg("aperture_y") = DipEdge::DEFAULT_aperture_y,
            py::arg("name") = py::none(),
            "Edge focusing associated with bend entry or exit."
        )
        .def_property("psi",
            [](DipEdge & dip_edge) { return dip_edge.m_psi; },
            [](DipEdge & dip_edge, amrex::ParticleReal psi) { dip_edge.m_psi = psi; },
            "Pole face angle in rad"
        )
        .def_property("rc",
            [](DipEdge & dip_edge) { return dip_edge.m_rc; },
            [](DipEdge & dip_edge, amrex::ParticleReal rc) { dip_edge.m_rc = rc; },
            "Radius of curvature in m"
        )
        .def_property("g",
            [](DipEdge & dip_edge) { return dip_edge.m_g; },
            [](DipEdge & dip_edge, amrex::ParticleReal g) { dip_edge.m_g = g; },
            "Gap parameter in m"
        )
        .def_property("R",
            [](DipEdge & dip_edge) { return dip_edge.m_R; },
            [](DipEdge & dip_edge, amrex::ParticleReal R) {
                if (R <= 0.0)
                     throw std::invalid_argument(R"(DipEdge parameter R must be > 0.)");
                dip_edge.m_R = R;
            },
            "Length scale for field integrals in m"
        )
        .def_property("K0",
            [](DipEdge & dip_edge) { return dip_edge.m_K0; },
            [](DipEdge & dip_edge, amrex::ParticleReal K0) { dip_edge.m_K0 = K0; },
            "Fringe field integral (unitless)"
        )
        .def_property("K1",
            [](DipEdge & dip_edge) { return dip_edge.m_K1; },
            [](DipEdge & dip_edge, amrex::ParticleReal K1) { dip_edge.m_K1 = K1; },
            "Fringe field integral (unitless)"
        )
        .def_property("K2",
            [](DipEdge & dip_edge) { return dip_edge.m_K2; },
            [](DipEdge & dip_edge, amrex::ParticleReal K2) { dip_edge.m_K2 = K2; },
            "Fringe field integral (unitless)"
        )
        .def_property("K3",
            [](DipEdge & dip_edge) { return dip_edge.m_K3; },
            [](DipEdge & dip_edge, amrex::ParticleReal K3) { dip_edge.m_K3 = K3; },
            "Fringe field integral (unitless)"
        )
        .def_property("K4",
            [](DipEdge & dip_edge) { return dip_edge.m_K4; },
            [](DipEdge & dip_edge, amrex::ParticleReal K4) { dip_edge.m_K4 = K4; },
            "Fringe field integral (unitless)"
        )
        .def_property("K5",
            [](DipEdge & dip_edge) { return dip_edge.m_K5; },
            [](DipEdge & dip_edge, amrex::ParticleReal K5) { dip_edge.m_K5 = K5; },
            "Fringe field integral (unitless)"
        )
        .def_property("K6",
            [](DipEdge & dip_edge) { return dip_edge.m_K6; },
            [](DipEdge & dip_edge, amrex::ParticleReal K6) { dip_edge.m_K6 = K6; },
            "Fringe field integral (unitless)"
        )
        .def_property("model",
            [](DipEdge & dip_edge) { return amrex::getEnumNameString(dip_edge.m_model); },
            [](DipEdge & dip_edge, std::string const & model) {
                dip_edge.m_model = amrex::getEnum<DipEdge::Model>(model);
            },
            "Fringe field model (linear or nonlinear)"
        )
        .def_property("location",
            [](DipEdge & dip_edge) { return amrex::getEnumNameString(dip_edge.m_location); },
            [](DipEdge & dip_edge, std::string const & location) {
                dip_edge.m_location = amrex::getEnum<DipEdge::Location>(location);
            },
            "Fringe field location (entry or exit)"
        )
        .def_property("modify_ref_part",
            [](DipEdge & dip_edge) { return dip_edge.m_modify_ref_part; },
            [](DipEdge & dip_edge, bool modify_ref_part) { dip_edge.m_modify_ref_part = modify_ref_part; },
            "Apply DipEdge to reference particle (boolean)."
        )

    ;
    register_push(py_DipEdge);
    register_reverse(py_DipEdge);

    py::class_<QuadEdge, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_QuadEdge(me, "QuadEdge");
    py_QuadEdge
        .def("__repr__",
             [](QuadEdge const & quadedge) {
                 return element_name(
                     quadedge,
                     std::make_pair("k", quadedge.m_k)
                 );
             }
        )
        .def("to_dict",
            [](QuadEdge const & quadedge) {
                return element_dict(
                    quadedge,
                    std::make_pair("k", quadedge.m_k),
                    std::make_pair("unit", quadedge.m_unit),
                    std::make_pair("flag", amrex::getEnumNameString(quadedge.m_flag))
                );
            }
        )
        .def(py::init([](
                amrex::ParticleReal k,
                int unit,
                std::string const & flag,
                amrex::ParticleReal dx,
                amrex::ParticleReal dy,
                amrex::ParticleReal rotation_degree,
                amrex::ParticleReal aperture_x,
                amrex::ParticleReal aperture_y,
                std::optional<std::string> name
             )
             {
                 QuadEdge::Location const fl = amrex::getEnum<QuadEdge::Location>(flag);
                 return new QuadEdge(k, unit, fl, dx, dy, rotation_degree, aperture_x, aperture_y, name);
             }),
             py::arg("k"),
             py::arg("unit") = QuadEdge::DEFAULT_unit,
             py::arg("flag") = amrex::getEnumNameString(QuadEdge::DEFAULT_flag),
             py::arg("dx") = QuadEdge::DEFAULT_dx,
             py::arg("dy") = QuadEdge::DEFAULT_dy,
             py::arg("rotation") = QuadEdge::DEFAULT_rotation_degree,
             py::arg("aperture_x") = QuadEdge::DEFAULT_aperture_x,
             py::arg("aperture_y") = QuadEdge::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             R"(A thin quadrupole fringe field element. Flag must be "entry" or "exit".)"
        )
        .def_property("k",
            [](QuadEdge & quadedge) { return quadedge.m_k; },
            [](QuadEdge & quadedge, amrex::ParticleReal k) { quadedge.m_k = k; },
            "quadrupole focusing strength (1/meter^2 OR T/m)"
        )
        .def_property("unit",
            [](QuadEdge & quadedge) { return quadedge.m_unit; },
            [](QuadEdge & quadedge, int unit) { quadedge.m_unit = unit; },
            "unit specification for quad strength"
        )
    ;
    register_push(py_QuadEdge);
    register_reverse(py_QuadEdge);

    py::class_<Drift, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_Drift(me, "Drift");
    py_Drift
        .def("__repr__",
             [](Drift const & drift) {
                 return element_name(drift);
             }
        )
        .def("to_dict",
             [](Drift const & drift) {
                 return element_dict(drift);
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("dx") = Drift::DEFAULT_dx,
             py::arg("dy") = Drift::DEFAULT_dy,
             py::arg("rotation") = Drift::DEFAULT_rotation_degree,
             py::arg("aperture_x") = Drift::DEFAULT_aperture_x,
             py::arg("aperture_y") = Drift::DEFAULT_aperture_y,
             py::arg("nslice") = Drift::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A drift."
        )
    ;
    register_push(py_Drift);
    register_reverse(py_Drift);

    py::class_<ExactDrift, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ExactDrift(me, "ExactDrift");
    py_ExactDrift
        .def("__repr__",
             [](ExactDrift const & exact_drift) {
                 return element_name(exact_drift);
             }
        )
        .def("to_dict",
             [](ExactDrift const & exact_drift) {
                 return element_dict(exact_drift);
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("dx") = ExactDrift::DEFAULT_dx,
             py::arg("dy") = ExactDrift::DEFAULT_dy,
             py::arg("rotation") = ExactDrift::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ExactDrift::DEFAULT_aperture_x,
             py::arg("aperture_y") = ExactDrift::DEFAULT_aperture_y,
             py::arg("nslice") = ExactDrift::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A Drift using the exact nonlinear map."
        )
    ;
    register_push(py_ExactDrift);
    register_reverse(py_ExactDrift);

    py::class_<ExactMultipole, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ExactMultipole(me, "ExactMultipole");
    py_ExactMultipole
        .def("__repr__",
             [](ExactMultipole const & exact_multipole) {
                 return element_name(
                     exact_multipole,
                     std::make_pair("unit", exact_multipole.m_unit)
                 );
             }
        )
        .def("to_dict",
             [](ExactMultipole const & exact_multipole) {
                 return element_dict(
                     exact_multipole,
                     std::make_pair("unit", exact_multipole.m_unit),
                     std::make_pair("k_normal", ExactMultipole::DynamicData::get(exact_multipole.m_id)->k_normal.host_const()),
                     std::make_pair("k_skew", ExactMultipole::DynamicData::get(exact_multipole.m_id)->k_skew.host_const()),
                     std::make_pair("mapsteps", exact_multipole.m_mapsteps),
                     std::make_pair("int_order", exact_multipole.m_int_order)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                std::vector<amrex::ParticleReal>,
                std::vector<amrex::ParticleReal>,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                int,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("k_normal"),
             py::arg("k_skew"),
             py::arg("unit") = ExactMultipole::DEFAULT_unit,
             py::arg("dx") = ExactMultipole::DEFAULT_dx,
             py::arg("dy") = ExactMultipole::DEFAULT_dy,
             py::arg("rotation") = ExactMultipole::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ExactMultipole::DEFAULT_aperture_x,
             py::arg("aperture_y") = ExactMultipole::DEFAULT_aperture_y,
             py::arg("int_order") = ExactMultipole::DEFAULT_int_order,
             py::arg("mapsteps") = ExactMultipole::DEFAULT_mapsteps,
             py::arg("nslice") = ExactMultipole::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A thick Multipole magnet using the exact nonlinear Hamiltonian."
        )
        .def_property("unit",
            [](ExactMultipole & exact_multipole) { return exact_multipole.m_unit; },
            [](ExactMultipole & exact_multipole, int unit) { exact_multipole.m_unit = unit; },
            "unit specification for multipole strength"
        )
        .def_property("int_order",
            [](ExactMultipole & exact_multipole) { return exact_multipole.m_int_order; },
            [](ExactMultipole & exact_multipole, int int_order) {
                if (int_order != 2 && int_order != 4 && int_order != 6)
                    throw std::invalid_argument("ExactMultipole: The order used for symplectic integration must be 2, 4 or 6.");
                exact_multipole.m_int_order = int_order;
            },
            "order of symplectic integration used for particle push in applied fields"
        )
        .def_property("mapsteps",
            [](ExactMultipole & exact_multipole) { return exact_multipole.m_mapsteps; },
            [](ExactMultipole & exact_multipole, int mapsteps) {
                if (mapsteps <= 0)
                    throw std::invalid_argument("ExactMultipole: mapsteps must be > 0!");
                exact_multipole.m_mapsteps = mapsteps;
            },
            "number of integration steps per slice used for particle push in the applied fields"
        )
    ;
    register_push(py_ExactMultipole);
    register_reverse(py_ExactMultipole);

    py::class_<ExactCFbend, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ExactCFbend(me, "ExactCFbend");
    py_ExactCFbend
        .def("__repr__",
             [](ExactCFbend const & exact_cfbend) {
                 return element_name(
                     exact_cfbend,
                     std::make_pair("unit", exact_cfbend.m_unit)
                 );
             }
        )
        .def("to_dict",
             [](ExactCFbend const & exact_cfbend) {
                 return element_dict(
                     exact_cfbend,
                     std::make_pair("unit", exact_cfbend.m_unit),
                     std::make_pair("k_normal", ExactCFbend::DynamicData::get(exact_cfbend.m_id)->k_normal.host_const()),
                     std::make_pair("k_skew", ExactCFbend::DynamicData::get(exact_cfbend.m_id)->k_skew.host_const()),
                     std::make_pair("mapsteps", exact_cfbend.m_mapsteps),
                     std::make_pair("int_order", exact_cfbend.m_int_order)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                std::vector<amrex::ParticleReal>,
                std::vector<amrex::ParticleReal>,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                int,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("k_normal"),
             py::arg("k_skew"),
             py::arg("unit") = ExactCFbend::DEFAULT_unit,
             py::arg("dx") = ExactCFbend::DEFAULT_dx,
             py::arg("dy") = ExactCFbend::DEFAULT_dy,
             py::arg("rotation") = ExactCFbend::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ExactCFbend::DEFAULT_aperture_x,
             py::arg("aperture_y") = ExactCFbend::DEFAULT_aperture_y,
             py::arg("int_order") = ExactCFbend::DEFAULT_int_order,
             py::arg("mapsteps") = ExactCFbend::DEFAULT_mapsteps,
             py::arg("nslice") = ExactCFbend::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A thick combined function bending magnet using the exact nonlinear Hamiltonian."
        )
        .def_property("unit",
            [](ExactCFbend & exact_cfbend) { return exact_cfbend.m_unit; },
            [](ExactCFbend & exact_cfbend, int unit) { exact_cfbend.m_unit = unit; },
            "unit specification for multipole strength"
        )
        .def_property("int_order",
            [](ExactCFbend & exact_cfbend) { return exact_cfbend.m_int_order; },
            [](ExactCFbend & exact_cfbend, int int_order) {
                if (int_order != 2 && int_order != 4 && int_order != 6)
                    throw std::invalid_argument("ExactCFbend: The order used for symplectic integration must be 2, 4 or 6.");
                exact_cfbend.m_int_order = int_order;
            },
            "order of symplectic integration used for particle push in applied fields"
        )
        .def_property("mapsteps",
            [](ExactCFbend & exact_cfbend) { return exact_cfbend.m_mapsteps; },
            [](ExactCFbend & exact_cfbend, int mapsteps) {
                if (mapsteps <= 0)
                    throw std::invalid_argument("ExactCFbend: mapsteps must be > 0!");
                exact_cfbend.m_mapsteps = mapsteps;
            },
            "number of integration steps per slice used for particle push in the applied fields"
        )
    ;
    register_push(py_ExactCFbend);
    register_reverse(py_ExactCFbend);

    py::class_<ExactQuad, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ExactQuad(me, "ExactQuad");
    py_ExactQuad
        .def("__repr__",
             [](ExactQuad const & exact_quad) {
                 return element_name(
                     exact_quad,
                     std::make_pair("k", exact_quad.m_k),
                     std::make_pair("unit", exact_quad.m_unit)
                 );
             }
        )
        .def("to_dict",
             [](ExactQuad const & exact_quad) {
                 return element_dict(
                     exact_quad,
                     std::make_pair("k", exact_quad.m_k),
                     std::make_pair("unit", exact_quad.m_unit),
                     std::make_pair("mapsteps", exact_quad.m_mapsteps),
                     std::make_pair("int_order", exact_quad.m_int_order)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                int,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("k"),
             py::arg("unit") = ExactQuad::DEFAULT_unit,
             py::arg("dx") = ExactQuad::DEFAULT_dx,
             py::arg("dy") = ExactQuad::DEFAULT_dy,
             py::arg("rotation") = ExactQuad::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ExactQuad::DEFAULT_aperture_x,
             py::arg("aperture_y") = ExactQuad::DEFAULT_aperture_y,
             py::arg("int_order") = ExactQuad::DEFAULT_int_order,
             py::arg("mapsteps") = ExactQuad::DEFAULT_mapsteps,
             py::arg("nslice") = ExactQuad::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A Quadrupole magnet using the exact nonlinear Hamiltonian."
        )
        .def_property("k",
            [](ExactQuad & exact_quad) { return exact_quad.m_k; },
            [](ExactQuad & exact_quad, amrex::ParticleReal k) { exact_quad.m_k = k; },
            "quadrupole strength in 1/m^2 (or T/m)"
        )
        .def_property("unit",
            [](ExactQuad & exact_quad) { return exact_quad.m_unit; },
            [](ExactQuad & exact_quad, int unit) { exact_quad.m_unit = unit; },
            "unit specification for quad strength"
        )
        .def_property("int_order",
            [](ExactQuad & exact_quad) { return exact_quad.m_int_order; },
            [](ExactQuad & exact_quad, int int_order) {
                if (int_order != 2 && int_order != 4 && int_order != 6)
                    throw std::invalid_argument("ExactQuad: The order used for symplectic integration must be 2, 4 or 6.");
                exact_quad.m_int_order = int_order;
            },
            "order of symplectic integration used for particle push in applied fields"
        )
        .def_property("mapsteps",
            [](ExactQuad & exact_quad) { return exact_quad.m_mapsteps; },
            [](ExactQuad & exact_quad, int mapsteps) {
                if (mapsteps <= 0)
                    throw std::invalid_argument("ExactQuad: mapsteps must be > 0!");
                exact_quad.m_mapsteps = mapsteps;
            },
            "number of integration steps per slice used for particle push in the applied fields"
        )
    ;
    register_push(py_ExactQuad);
    register_reverse(py_ExactQuad);

    py::class_<ExactSbend, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ExactSbend(me, "ExactSbend");
    py_ExactSbend
        .def("__repr__",
             [](ExactSbend const & exact_sbend) {
                 return element_name(
                     exact_sbend,
                     std::make_pair("phi", exact_sbend.m_phi / ExactSbend::degree2rad),
                     std::make_pair("B", exact_sbend.m_B)
                 );
             }
        )
        .def("to_dict",
            [](ExactSbend const & exact_sbend, bool in_degrees) {
                if (in_degrees) {
                    return element_dict(
                        exact_sbend,
                        std::make_pair("phi", exact_sbend.m_phi / ExactSbend::degree2rad),
                                                                   // once fixed, update src/python/impactx/extensions/KnownElementsList.py
                        std::make_pair("B", exact_sbend.m_B)
                    );
                } else {
                    // legacy: buggy radians instead of degrees
                    py::warnings::warn(
                        "Warning: ExactSbend.to_dict() has a known bug, "
                        "returning phi in radians than degrees. "
                        "Please use ExactSbend.to_dict(in_degrees=True) instead.",
                        PyExc_RuntimeWarning,
                        2
                    );
                    return element_dict(
                        exact_sbend,
                        std::make_pair("phi", exact_sbend.m_phi),  // BUG: constructor is in degrees
                        std::make_pair("B", exact_sbend.m_B)
                    );
                }
            },
            py::arg("in_degrees") = false
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("phi"),
             py::arg("B") = ExactSbend::DEFAULT_B,
             py::arg("dx") = ExactSbend::DEFAULT_dx,
             py::arg("dy") = ExactSbend::DEFAULT_dy,
             py::arg("rotation") = ExactSbend::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ExactSbend::DEFAULT_aperture_x,
             py::arg("aperture_y") = ExactSbend::DEFAULT_aperture_y,
             py::arg("nslice") = ExactSbend::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "An ideal sector bend using the exact nonlinear map.  When B = 0, the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B."
        )
        .def("rc", &ExactSbend::rc,
            py::arg("ref"),
            "Radius of curvature in m"
        )
        .def_property("phi",
            [](ExactSbend & exact_sbend) { return exact_sbend.m_phi; },
            [](ExactSbend & exact_sbend, amrex::ParticleReal phi) { exact_sbend.m_phi = phi; },
            "Bend angle in radian"
        )
        /* BUG, should be in degrees like this:
        .def_property("phi",
            [](ExactSbend & exact_sbend) { return exact_sbend.m_phi / ExactSbend::degree2rad; },
            [](ExactSbend & exact_sbend, amrex::ParticleReal phi_deg) {
                exact_sbend.m_phi = phi_deg * ExactSbend::degree2rad;
            },
            "Bend angle in degrees"
        )
        */
        .def_property("B",
            [](ExactSbend & exact_sbend) { return exact_sbend.m_B; },
            [](ExactSbend & exact_sbend, amrex::ParticleReal B) { exact_sbend.m_B = B; },
            "Magnetic field in Tesla; when B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B"
        )
    ;
    register_push(py_ExactSbend);
    register_reverse(py_ExactSbend);

    py::class_<Kicker, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_Kicker(me, "Kicker");
    py_Kicker
        .def("__repr__",
             [](Kicker const & kicker) {
                 return element_name(
                     kicker,
                     std::make_pair("xkick", kicker.m_xkick),
                     std::make_pair("ykick", kicker.m_ykick)
                 );
             }
        )
        .def("to_dict",
            [](Kicker const & kicker) {
                return element_dict(
                    kicker,
                    std::make_pair("xkick", kicker.m_xkick),
                    std::make_pair("ykick", kicker.m_ykick),
                    std::make_pair("unit", Kicker::unit_name(kicker.m_unit))
                );
            }
        )
        .def(py::init([](
                amrex::ParticleReal xkick,
                amrex::ParticleReal ykick,
                std::string const & unit,
                amrex::ParticleReal dx,
                amrex::ParticleReal dy,
                amrex::ParticleReal rotation_degree,
                amrex::ParticleReal aperture_x,
                amrex::ParticleReal aperture_y,
                std::optional<std::string> name
             )
             {
                 Kicker::UnitSystem const u = Kicker::unit_from_name(unit);
                 return new Kicker(xkick, ykick, u, dx, dy, rotation_degree, aperture_x, aperture_y, name);
             }),
             py::arg("xkick"),
             py::arg("ykick"),
             py::arg("unit") = Kicker::unit_name(Kicker::DEFAULT_unit),
             py::arg("dx") = Kicker::DEFAULT_dx,
             py::arg("dy") = Kicker::DEFAULT_dy,
             py::arg("rotation") = Kicker::DEFAULT_rotation_degree,
             py::arg("aperture_x") = Kicker::DEFAULT_aperture_x,
             py::arg("aperture_y") = Kicker::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             R"(A thin transverse kicker element. Kicks are for unit "dimensionless" or in "T-m".)"
        )
        .def_property("xkick",
            [](Kicker & kicker) { return kicker.m_xkick; },
            [](Kicker & kicker, amrex::ParticleReal xkick) { kicker.m_xkick = xkick; },
            "horizontal kick strength (dimensionless OR T-m)"
        )
        .def_property("ykick",
            [](Kicker & kicker) { return kicker.m_ykick; },
            [](Kicker & kicker, amrex::ParticleReal ykick) { kicker.m_ykick = ykick; },
            "vertical kick strength (dimensionless OR T-m)"
        )
        // TODO unit
    ;
    register_push(py_Kicker);
    register_reverse(py_Kicker);

    py::class_<Multipole, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_Multipole(me, "Multipole");
    py_Multipole
        .def("__repr__",
             [](Multipole const & multipole) {
                 return element_name(
                     multipole,
                     std::make_pair("multipole", multipole.m_multipole),
                     std::make_pair("K_normal", multipole.m_Kn),
                     std::make_pair("K_skew", multipole.m_Ks)
                 );
             }
        )
        .def("to_dict",
            [](Multipole const & multipole) {
                return element_dict(
                    multipole,
                    std::make_pair("multipole", multipole.m_multipole),
                    std::make_pair("K_normal", multipole.m_Kn),
                    std::make_pair("K_skew", multipole.m_Ks)
                );
            }
        )
        .def(py::init<
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("multipole"),
             py::arg("K_normal"),
             py::arg("K_skew"),
             py::arg("dx") = Multipole::DEFAULT_dx,
             py::arg("dy") = Multipole::DEFAULT_dy,
             py::arg("rotation") = Multipole::DEFAULT_rotation_degree,
             py::arg("aperture_x") = Multipole::DEFAULT_aperture_x,
             py::arg("aperture_y") = Multipole::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             "A general thin multipole element."
        )
        .def_property("multipole",
            [](Multipole & multipole) { return multipole.m_multipole; },
            [](Multipole & multipole, int multipole_index) {
                multipole.m_multipole = multipole_index;
                multipole.compute_factorial();
            },
            "index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)"
        )
        .def_property("K_normal",
            [](Multipole & multipole) { return multipole.m_Kn; },
            [](Multipole & multipole, amrex::ParticleReal K_normal) { multipole.m_Kn = K_normal; },
            "Integrated normal multipole coefficient (1/meter^m)"
        )
        .def_property("K_skew",
            [](Multipole & multipole) { return multipole.m_Ks; },
            [](Multipole & multipole, amrex::ParticleReal K_skew) { multipole.m_Ks = K_skew; },
            "Integrated skew multipole coefficient (1/meter^m)"
        )
    ;
    register_push(py_Multipole);
    register_reverse(py_Multipole);

    py::class_<Empty, elements::mixin::Named, elements::mixin::Thin> py_Empty(me, "Empty");
    py_Empty
        .def("__repr__",
             [](Empty const & /* empty */) {
                 return std::string("<impactx.elements.Empty>");
             }
        )
        .def("to_dict",
            [](Empty const & empty) {
                return element_dict(empty);
            }
        )
        .def(py::init<>(),
             "This element does nothing."
        )
    ;
    register_push(py_Empty);
    register_reverse(py_Empty);

    py::class_<Marker, elements::mixin::Named, elements::mixin::Thin> py_Marker(me, "Marker");
    py_Marker
        .def("__repr__",
             [](Marker const & marker) {
                 return element_name(marker);
             }
        )
        .def("to_dict",
            [](Marker const & marker) {
                return element_dict(marker);
            }
        )
        .def(py::init<std::string>(),
             py::arg("name"),
             "This named element does nothing."
        )
    ;
    register_push(py_Marker);
    register_reverse(py_Marker);

    py::class_<NonlinearLens, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_NonlinearLens(me, "NonlinearLens");
    py_NonlinearLens
        .def("__repr__",
             [](NonlinearLens const & nl) {
                 return element_name(
                     nl,
                     std::make_pair("knll", nl.m_knll),
                     std::make_pair("cnll", nl.m_cnll)
                 );
             }
        )
        .def("to_dict",
            [](NonlinearLens const & nl) {
                return element_dict(
                    nl,
                    std::make_pair("knll", nl.m_knll),
                    std::make_pair("cnll", nl.m_cnll)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("knll"),
             py::arg("cnll"),
             py::arg("dx") = NonlinearLens::DEFAULT_dx,
             py::arg("dy") = NonlinearLens::DEFAULT_dy,
             py::arg("rotation") = NonlinearLens::DEFAULT_rotation_degree,
             py::arg("aperture_x") = NonlinearLens::DEFAULT_aperture_x,
             py::arg("aperture_y") = NonlinearLens::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             "Single short segment of the nonlinear magnetic insert element."
        )
        .def_property("knll",
            [](NonlinearLens & nl) { return nl.m_knll; },
            [](NonlinearLens & nl, amrex::ParticleReal knll) { nl.m_knll = knll; },
            "integrated strength of the nonlinear lens (m)"
        )
        .def_property("cnll",
            [](NonlinearLens & nl) { return nl.m_cnll; },
            [](NonlinearLens & nl, amrex::ParticleReal cnll) { nl.m_cnll = cnll; },
            "distance of singularities from the origin (m)"
        )
    ;
    register_push(py_NonlinearLens);
    register_reverse(py_NonlinearLens);

    py::class_<PlaneXYRot, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment> py_PlaneXYRot(me, "PlaneXYRot");
    py_PlaneXYRot
        .def("__repr__",
             [](PlaneXYRot const & plane_xyrot) {
                 return element_name(
                     plane_xyrot,
                     std::make_pair("angle", plane_xyrot.m_phi / PlaneXYRot::degree2rad)
                 );
             }
        )
        .def("to_dict",
            [](PlaneXYRot const & plane_xyrot, bool in_degrees) {
                if (in_degrees) {
                    return element_dict(
                        plane_xyrot,
                        std::make_pair("angle", plane_xyrot.m_phi / PlaneXYRot::degree2rad)
                                                                      // once fixed, update src/python/impactx/extensions/KnownElementsList.py
                    );
                } else {
                    // legacy: buggy radians instead of degrees
                    py::warnings::warn(
                        "Warning: PlaneXYRot.to_dict() has a known bug, "
                        "returning angle in radians than degrees. "
                        "Please use PlaneXYRot.to_dict(in_degrees=True) instead.",
                        PyExc_RuntimeWarning,
                        2
                    );
                    return element_dict(
                        plane_xyrot,
                        std::make_pair("angle", plane_xyrot.m_phi)  // BUG: constructor is in degrees
                    );
                }
            },
            py::arg("in_degrees") = false
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("angle"),
             py::arg("dx") = PlaneXYRot::DEFAULT_dx,
             py::arg("dy") = PlaneXYRot::DEFAULT_dy,
             py::arg("rotation") = PlaneXYRot::DEFAULT_rotation_degree,
             py::arg("name") = py::none(),
             "A rotation in the x-y plane."
        )
        .def_property("angle",
            [](PlaneXYRot & plane_xyrot) { return plane_xyrot.m_phi; },
            [](PlaneXYRot & plane_xyrot, amrex::ParticleReal phi) { plane_xyrot.m_phi = phi; },
            "Rotation angle (rad)."
        )
        /* BUG: should be in degrees
        .def_property("angle",
            [](PlaneXYRot & plane_xyrot) { return plane_xyrot.m_phi / PlaneXYRot::degree2rad; },
            [](PlaneXYRot & plane_xyrot, amrex::ParticleReal phi_deg) {
                plane_xyrot.m_phi = phi_deg * PlaneXYRot::degree2rad;
            },
            "Rotation angle in degrees (in the x-y plane about the reference trajectory)."
        )
        */
    ;
    register_push(py_PlaneXYRot);
    register_reverse(py_PlaneXYRot);

    py::class_<PolygonAperture, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment> py_PolygonAperture(me, "PolygonAperture");
    py_PolygonAperture
        .def("__repr__",
             [](PolygonAperture const & polygon_aperture) {
                 return element_name(
                    polygon_aperture,
                    std::make_pair("action", polygon_aperture.action_name(polygon_aperture.m_action))
                );
             }
        )
        .def("to_dict",
            [](PolygonAperture const & polygon_aperture) {
                using namespace amrex::literals;
                return element_dict(
                    polygon_aperture,
                    std::make_pair("vertices_x", PolygonAperture::DynamicData::get(polygon_aperture.m_id)->x.host_const()),
                    std::make_pair("vertices_y", PolygonAperture::DynamicData::get(polygon_aperture.m_id)->y.host_const()),
                    std::make_pair("min_radius2", polygon_aperture.m_min_radius2),
                    std::make_pair("action", polygon_aperture.action_name(polygon_aperture.m_action)),
                    std::make_pair("repeat_x", polygon_aperture.m_repeat_x),
                    std::make_pair("repeat_y", polygon_aperture.m_repeat_y),
                    std::make_pair("shift_odd_x", polygon_aperture.m_shift_odd_x)
                );
            }
        )
        .def(py::init([](
                 std::vector<amrex::ParticleReal> vertices_x,
                 std::vector<amrex::ParticleReal> vertices_y,
                 amrex::ParticleReal min_radius2,
                 amrex::ParticleReal repeat_x,
                 amrex::ParticleReal repeat_y,
                 bool shift_odd_x,
                 std::string const & action,
                 amrex::ParticleReal dx,
                 amrex::ParticleReal dy,
                 amrex::ParticleReal rotation_degree,
                 std::optional<std::string> name
             )
             {
                 return new PolygonAperture(
                     vertices_x, vertices_y, min_radius2, repeat_x, repeat_y, shift_odd_x,
                     amrex::getEnum<PolygonAperture::Action>(action),
                     dx, dy, rotation_degree, name
                 );
             }),
             py::arg("vertices_x"),
             py::arg("vertices_y"),
             py::arg("min_radius2") = PolygonAperture::DEFAULT_min_radius2,
             py::arg("repeat_x") = PolygonAperture::DEFAULT_repeat_x,
             py::arg("repeat_y") = PolygonAperture::DEFAULT_repeat_y,
             py::arg("shift_odd_x") = PolygonAperture::DEFAULT_shift_odd_x,
             py::arg("action") = amrex::getEnumNameString(PolygonAperture::DEFAULT_action),
             py::arg("dx") = PolygonAperture::DEFAULT_dx,
             py::arg("dy") = PolygonAperture::DEFAULT_dy,
             py::arg("rotation") = PolygonAperture::DEFAULT_rotation_degree,
             py::arg("name") = py::none(),
             "A short collimator element described by a polygon with vertices given by their x and y coordinates."
        )
        .def_property("action",
            [](PolygonAperture & ap)
            {
                return ap.action_name(ap.m_action);
            },
            [](PolygonAperture & ap, std::string const & action)
            {
                ap.m_action = amrex::getEnum<PolygonAperture::Action>(action);
            },
            "action type (transmit, absorb)"
        )
        .def_property("min_radius2",
            [](PolygonAperture & pa) { return pa.m_min_radius2; },
            [](PolygonAperture & pa, amrex::ParticleReal mr2) { pa.m_min_radius2 = mr2; },
            "All particles with radius squared smaller than min_radius2 pass the aperture"
        )
        .def_property("repeat_x",
            [](PolygonAperture & ap) { return ap.m_repeat_x; },
            [](PolygonAperture & ap, amrex::ParticleReal repeat_x) { ap.m_repeat_x = repeat_x; },
            "horizontal period for repeated aperture masking"
        )
        .def_property("repeat_y",
            [](PolygonAperture & ap) { return ap.m_repeat_y; },
            [](PolygonAperture & ap, amrex::ParticleReal repeat_y) { ap.m_repeat_y = repeat_y; },
            "vertical period for repeated aperture masking"
        )
        .def_property("shift_odd_x",
            [](PolygonAperture & ap) { return ap.m_shift_odd_x; },
            [](PolygonAperture & ap, bool shift_odd_x) { ap.m_shift_odd_x = shift_odd_x; },
            "for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. "
            "Use alignment offsets dx,dy to move whole mask as needed."
        )
    ;
    register_push(py_PolygonAperture);
    register_reverse(py_PolygonAperture);

    py::class_<Programmable, elements::mixin::Named>(me, "Programmable", py::dynamic_attr())
        .def("__repr__",
             [](Programmable const & prg) {
                 return element_name(prg);
             }
        )
        .def("to_dict",
            [](Programmable const & prg) {
                return element_dict(
                    prg,
                    std::make_pair("ds", prg.m_ds),
                    std::make_pair("nslice", prg.m_nslice)
                );
            }
        )
        .def(py::init<
                 amrex::ParticleReal,
                 int,
                 std::optional<std::string>
             >(),
             py::arg("ds") = Programmable::DEFAULT_ds,
             py::arg("nslice") = Programmable::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A programmable beam optics element."
        )
        .def_property("nslice",
            [](Programmable & p) { return p.nslice(); },
            [](Programmable & p, int nslice) {
                if (nslice <= 0)
                    throw std::invalid_argument("Programmable: nslice must be > 0!");
                p.m_nslice = nslice;
            }
        )
        .def_property("ds",
              [](Programmable & p) { return p.ds(); },
              [](Programmable & p, amrex::ParticleReal ds) { p.m_ds = ds; }
        )
        .def_property("threadsafe",
            [](Programmable & p) { return p.m_threadsafe; },
            [](Programmable & p, bool threadsafe) { p.m_threadsafe = threadsafe; },
            "allow threading via OpenMP for the particle iterator loop, default=False (note: if OMP backend is active)"
        )
        .def_property("push",
              [](Programmable & p) { return p.m_push; },
              [](Programmable & p,
                 std::function<void(ImpactXParticleContainer *, int, int)> new_hook
              ) { p.m_push = std::move(new_hook); },
              "hook for push of whole container (pc, step, period)"
        )
        .def_property("beam_particles",
              [](Programmable & p) { return p.m_beam_particles; },
              [](Programmable & p,
                 std::function<void(ImpactXParticleContainer::iterator *, RefPart &)> new_hook
              ) { p.m_beam_particles = std::move(new_hook); },
              "hook for beam particles (pti, RefPart)"
        )
        .def_property("ref_particle",
              [](Programmable & p) { return p.m_ref_particle; },
              [](Programmable & p,
                 std::function<void(RefPart &)> new_hook
              ) { p.m_ref_particle = std::move(new_hook); },
              "hook for reference particle (RefPart)"
        )
    ;

    py::class_<Quad, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_Quad(me, "Quad");
    py_Quad
        .def("__repr__",
             [](Quad const & quad) {
                 return element_name(
                     quad,
                     std::make_pair("k", quad.m_k)
                 );
             }
        )
        .def("to_dict",
             [](Quad const & quad) {
                 return element_dict(
                     quad,
                     std::make_pair("k", quad.m_k)
                 );
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("k"),
             py::arg("dx") = Quad::DEFAULT_dx,
             py::arg("dy") = Quad::DEFAULT_dy,
             py::arg("rotation") = Quad::DEFAULT_rotation_degree,
             py::arg("aperture_x") = Quad::DEFAULT_aperture_x,
             py::arg("aperture_y") = Quad::DEFAULT_aperture_y,
             py::arg("nslice") = Quad::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A Quadrupole magnet."
        )
        .def_property("k",
            [](Quad & quad) { return quad.m_k; },
            [](Quad & quad, amrex::ParticleReal k) { quad.m_k = k; },
            "Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing"
        )
    ;
    register_push(py_Quad);
    register_reverse(py_Quad);

    py::class_<RFCavity, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_RFCavity(me, "RFCavity");
    py_RFCavity
        .def("__repr__",
             [](RFCavity const & rfc) {
                 return element_name(
                     rfc,
                     std::make_pair("escale", rfc.m_escale),
                     std::make_pair("freq", rfc.m_freq),
                     std::make_pair("phase", rfc.m_phase)
                 );
             }
        )
        .def("to_dict",
            [](RFCavity const & rfc) {
                return element_dict(
                    rfc,
                    std::make_pair("escale", rfc.m_escale),
                    std::make_pair("freq", rfc.m_freq),
                    std::make_pair("phase", rfc.m_phase),
                    std::make_pair("cos_coefficients", RFCavity::DynamicData::get(rfc.m_id)->cos.host_const()),
                    std::make_pair("sin_coefficients", RFCavity::DynamicData::get(rfc.m_id)->sin.host_const()),
                    std::make_pair("z_data", RFCavity::WakeData::get(rfc.m_wake_id)->z_data.host_const()),
                    std::make_pair("wake_x_data", RFCavity::WakeData::get(rfc.m_wake_id)->wake_x_data.host_const()),
                    std::make_pair("wake_y_data", RFCavity::WakeData::get(rfc.m_wake_id)->wake_y_data.host_const()),
                    std::make_pair("wake_z_data", RFCavity::WakeData::get(rfc.m_wake_id)->wake_z_data.host_const()),
                    std::make_pair("mapsteps", rfc.m_mapsteps)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::vector<amrex::ParticleReal>,
                std::vector<amrex::ParticleReal>,
                std::optional<std::vector<amrex::ParticleReal>>,
                std::optional<std::vector<amrex::ParticleReal>>,
                std::optional<std::vector<amrex::ParticleReal>>,
                std::optional<std::vector<amrex::ParticleReal>>,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("escale"),
             py::arg("freq"),
             py::arg("phase"),
             py::arg("cos_coefficients"),
             py::arg("sin_coefficients"),
             py::arg("z_data") = py::none(),
             py::arg("wake_x_data") = py::none(),
             py::arg("wake_y_data") = py::none(),
             py::arg("wake_z_data") = py::none(),
             py::arg("dx") = RFCavity::DEFAULT_dx,
             py::arg("dy") = RFCavity::DEFAULT_dy,
             py::arg("rotation") = RFCavity::DEFAULT_rotation_degree,
             py::arg("aperture_x") = RFCavity::DEFAULT_aperture_x,
             py::arg("aperture_y") = RFCavity::DEFAULT_aperture_y,
             py::arg("mapsteps") = RFCavity::DEFAULT_mapsteps,
             py::arg("nslice") = RFCavity::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "An RF cavity."
        )
        .def_property("escale",
            [](RFCavity & rfc) { return rfc.m_escale; },
            [](RFCavity & rfc, amrex::ParticleReal escale) { rfc.m_escale = escale; },
            "scaling factor for on-axis RF electric field in 1/m = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV)"
        )
        .def_property("freq",
            [](RFCavity & rfc) { return rfc.m_freq; },
            [](RFCavity & rfc, amrex::ParticleReal freq) { rfc.m_freq = freq; },
            "RF frequency in Hz"
        )
        .def_property("phase",
            [](RFCavity & rfc) { return rfc.m_phase; },
            [](RFCavity & rfc, amrex::ParticleReal phase) { rfc.m_phase = phase; },
            "RF driven phase in degrees"
        )
        // TODO cos_coefficients
        // TODO sin_coefficients
        .def_property("mapsteps",
            [](RFCavity & rfc) { return rfc.m_mapsteps; },
            [](RFCavity & rfc, int mapsteps) {
                if (mapsteps <= 0)
                    throw std::invalid_argument("RFCavity: mapsteps must be > 0!");
                rfc.m_mapsteps = mapsteps;
            },
            "number of integration steps per slice used for map and reference particle push in applied fields"
        )
        .def_property_readonly("map",
            [](RFCavity const & rfc) { return rfc.m_map; },
            "linearized transport map around the reference particle (valid after a reference-particle push)"
        )
        .def_property_readonly("spin_coupling",
            [](RFCavity const & rfc) { return rfc.m_spin_coupling; },
            "linearized spin-orbit coupling matrix (valid after a reference-particle push)"
        )
    ;
    register_push(py_RFCavity);
    register_reverse(py_RFCavity);

    py::class_<Sbend, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_Sbend(me, "Sbend");
    py_Sbend
        .def("__repr__",
             [](Sbend const & sbend) {
                 return element_name(
                     sbend,
                     std::make_pair("rc", sbend.m_rc)
                 );
             }
        )
        .def("to_dict",
            [](Sbend const & sbend) {
                return element_dict(
                    sbend,
                    std::make_pair("rc", sbend.m_rc)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("rc"),
             py::arg("dx") = Sbend::DEFAULT_dx,
             py::arg("dy") = Sbend::DEFAULT_dy,
             py::arg("rotation") = Sbend::DEFAULT_rotation_degree,
             py::arg("aperture_x") = Sbend::DEFAULT_aperture_x,
             py::arg("aperture_y") = Sbend::DEFAULT_aperture_y,
             py::arg("nslice") = Sbend::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "An ideal sector bend."
        )
        .def("rc", &Sbend::rc,
            py::arg("ref") = py::none(),
            "Radius of curvature in m"
        )
    ;
    register_push(py_Sbend);
    register_reverse(py_Sbend);

    py::class_<CFbend, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_CFbend(me, "CFbend");
    py_CFbend
        .def("__repr__",
             [](CFbend const & cfbend) {
                 return element_name(
                     cfbend,
                     std::make_pair("rc", cfbend.m_rc),
                     std::make_pair("k", cfbend.m_k)
                 );
             }
        )
        .def("to_dict",
            [](CFbend const & cfbend) {
                return element_dict(
                    cfbend,
                    std::make_pair("rc", cfbend.m_rc),
                    std::make_pair("k", cfbend.m_k)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("rc"),
             py::arg("k"),
             py::arg("dx") = CFbend::DEFAULT_dx,
             py::arg("dy") = CFbend::DEFAULT_dy,
             py::arg("rotation") = CFbend::DEFAULT_rotation_degree,
             py::arg("aperture_x") = CFbend::DEFAULT_aperture_x,
             py::arg("aperture_y") = CFbend::DEFAULT_aperture_y,
             py::arg("nslice") = CFbend::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "An ideal combined function bend (sector bend with quadrupole component)."
        )
        .def_property("rc",
            [](CFbend & cfbend) { return cfbend.m_rc; },
            [](CFbend & cfbend, amrex::ParticleReal rc) { cfbend.m_rc = rc; },
            "Radius of curvature in m"
        )
        .def_property("k",
            [](CFbend & cfbend) { return cfbend.m_k; },
            [](CFbend & cfbend, amrex::ParticleReal k) { cfbend.m_k = k; },
            "Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing"
        )
    ;
    register_push(py_CFbend);
    register_reverse(py_CFbend);

    py::class_<Buncher, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_Buncher(me, "Buncher");
    py_Buncher
        .def("__repr__",
             [](Buncher const & buncher) {
                 return element_name(
                     buncher,
                     std::make_pair("V", buncher.m_V),
                     std::make_pair("k", buncher.m_k)
                 );
             }
        )
        .def("to_dict",
            [](Buncher const & buncher) {
                return element_dict(
                    buncher,
                    std::make_pair("V", buncher.m_V),
                    std::make_pair("k", buncher.m_k)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("V"),
             py::arg("k"),
             py::arg("dx") = Buncher::DEFAULT_dx,
             py::arg("dy") = Buncher::DEFAULT_dy,
             py::arg("rotation") = Buncher::DEFAULT_rotation_degree,
             py::arg("aperture_x") = Buncher::DEFAULT_aperture_x,
             py::arg("aperture_y") = Buncher::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             "A short linear RF cavity element at zero-crossing for bunching."
        )
        .def_property("V",
            [](Buncher & buncher) { return buncher.m_V; },
            [](Buncher & buncher, amrex::ParticleReal V) { buncher.m_V = V; },
            "Normalized RF voltage drop V = Emax*L/(c*Brho)"
        )
        .def_property("k",
            [](Buncher & buncher) { return buncher.m_k; },
            [](Buncher & buncher, amrex::ParticleReal k) { buncher.m_k = k; },
            "Wavenumber of RF in 1/m"
        )
    ;
    register_push(py_Buncher);
    register_reverse(py_Buncher);

    py::class_<ShortRF, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ShortRF(me, "ShortRF");
    py_ShortRF
        .def("__repr__",
             [](ShortRF const & short_rf) {
                 return element_name(
                     short_rf,
                     std::make_pair("V", short_rf.m_V),
                     std::make_pair("freq", short_rf.m_freq),
                     std::make_pair("phase", short_rf.m_phase)
                 );
             }
        )
        .def("to_dict",
            [](ShortRF const & short_rf) {
                return element_dict(
                    short_rf,
                    std::make_pair("V", short_rf.m_V),
                    std::make_pair("freq", short_rf.m_freq),
                    std::make_pair("phase", short_rf.m_phase)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("V"),
             py::arg("freq"),
             py::arg("phase") = ShortRF::DEFAULT_phase,
             py::arg("dx") = ShortRF::DEFAULT_dx,
             py::arg("dy") = ShortRF::DEFAULT_dy,
             py::arg("rotation") = ShortRF::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ShortRF::DEFAULT_aperture_x,
             py::arg("aperture_y") = ShortRF::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             "A short RF cavity element."
        )
        .def_property("V",
            [](ShortRF & short_rf) { return short_rf.m_V; },
            [](ShortRF & short_rf, amrex::ParticleReal V) { short_rf.m_V = V; },
            "Normalized RF voltage V = maximum energy gain/(m*c^2)"
        )
        .def_property("freq",
            [](ShortRF & short_rf) { return short_rf.m_freq; },
            [](ShortRF & short_rf, amrex::ParticleReal freq) { short_rf.m_freq = freq; },
            "RF frequency in Hz"
        )
        .def_property("phase",
            [](ShortRF & short_rf) { return short_rf.m_phase; },
            [](ShortRF & short_rf, amrex::ParticleReal phase) { short_rf.m_phase = phase; },
            "RF synchronous phase in degrees (phase = 0 corresponds to maximum energy gain, phase = -90 corresponds go zero energy gain for bunching)"
        )
    ;
    register_push(py_ShortRF);
    register_reverse(py_ShortRF);

    py::class_<SoftSolenoid, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_SoftSolenoid(me, "SoftSolenoid");
    py_SoftSolenoid
        .def("__repr__",
             [](SoftSolenoid const & soft_sol) {
                 return element_name(
                     soft_sol,
                     std::make_pair("bscale", soft_sol.m_bscale)
                 );
             }
        )
        .def("to_dict",
            [](SoftSolenoid const & soft_sol) {
                return element_dict(
                    soft_sol,
                    std::make_pair("bscale", soft_sol.m_bscale),
                    std::make_pair("unit", soft_sol.m_unit),
                    std::make_pair("cos_coefficients", SoftSolenoid::DynamicData::get(soft_sol.m_id)->cos.host_const()),
                    std::make_pair("sin_coefficients", SoftSolenoid::DynamicData::get(soft_sol.m_id)->sin.host_const()),
                    std::make_pair("mapsteps", soft_sol.m_mapsteps)
                );
            }
        )
        .def(py::init<
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 std::vector<amrex::ParticleReal>,
                 std::vector<amrex::ParticleReal>,
                 int,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 int,
                 int,
                 std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("bscale"),
             py::arg("cos_coefficients"),
             py::arg("sin_coefficients"),
             py::arg("unit") = SoftSolenoid::DEFAULT_unit,
             py::arg("dx") = SoftSolenoid::DEFAULT_dx,
             py::arg("dy") = SoftSolenoid::DEFAULT_dy,
             py::arg("rotation") = SoftSolenoid::DEFAULT_rotation_degree,
             py::arg("aperture_x") = SoftSolenoid::DEFAULT_aperture_x,
             py::arg("aperture_y") = SoftSolenoid::DEFAULT_aperture_y,
             py::arg("mapsteps") = SoftSolenoid::DEFAULT_mapsteps,
             py::arg("nslice") = SoftSolenoid::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A soft-edge solenoid."
        )
        .def_property("bscale",
            [](SoftSolenoid & soft_sol) { return soft_sol.m_bscale; },
            [](SoftSolenoid & soft_sol, amrex::ParticleReal bscale) { soft_sol.m_bscale = bscale; },
            "Scaling factor for on-axis magnetic field Bz in inverse meters (if unit = 0) or magnetic field Bz in T (SI units, if unit = 1)"
        )
        /* TODO Before we can expose this we need to ensure that sim.lattice takes a list of shared pointers
         *    and not copies of elements. Otherwise, users can invalidate their cached host pointers with
         *      sol = SoftSolenoid(...)
         *      sim.lattice.append(sol)  # copy in lattice has same m_id
         *      sol.cos_coef = new_values  # updates data of m_id and sol.m_cos_h_data, but NOT the lattice copy
        .def_property("cos_coefficients",
            [](SoftSolenoid & soft_sol) {
                return SoftSolenoid::DynamicData::get(soft_sol.m_id)->h_cos;
            },
            [](SoftSolenoid & soft_sol, std::vector<amrex::ParticleReal> v) {
                auto & coef = *SoftSolenoid::DynamicData::get(soft_sol.m_id);
                coef.h_cos = std::move(v);
                coef.mark_dirty();
                soft_sol.m_cos_h_data = coef.h_cos.data();
            },
            "cosine coefficients in Fourier expansion of on-axis magnetic field Bz"
        )
        .def_property("sin_coefficients",
            [](SoftSolenoid & soft_sol) {
                return SoftSolenoid::DynamicData::get(soft_sol.m_id)->h_sin;
            },
            [](SoftSolenoid & soft_sol, std::vector<amrex::ParticleReal> v) {
                auto & coef = *SoftSolenoid::DynamicData::get(soft_sol.m_id);
                coef.h_sin = std::move(v);
                coef.mark_dirty();
                soft_sol.m_sin_h_data = coef.h_sin.data();
            },
            "sine coefficients in Fourier expansion of on-axis magnetic field Bz"
        )
        */
        .def_property("unit",
            [](SoftSolenoid & soft_sol) { return soft_sol.m_unit; },
            [](SoftSolenoid & soft_sol, int unit) { soft_sol.m_unit = unit; },
            "specification of units for scaling of the on-axis longitudinal magnetic field"
        )
        .def_property("mapsteps",
            [](SoftSolenoid & soft_sol) { return soft_sol.m_mapsteps; },
            [](SoftSolenoid & soft_sol, int mapsteps) {
                if (mapsteps <= 0)
                    throw std::invalid_argument("SoftSolenoid: mapsteps must be > 0!");
                soft_sol.m_mapsteps = mapsteps;
            },
            "number of integration steps per slice used for map and reference particle push in applied fields"
        )
        .def_property_readonly("map",
            [](SoftSolenoid const & soft_sol) { return soft_sol.m_map; },
            "linearized transport map around the reference particle (valid after a reference-particle push)"
        )
        .def_property_readonly("spin_coupling",
            [](SoftSolenoid const & soft_sol) { return soft_sol.m_spin_coupling; },
            "linearized spin-orbit coupling matrix (valid after a reference-particle push)"
        )
        .def_property_readonly("spin_rotation_vector",
            [](SoftSolenoid const & soft_sol) { return soft_sol.m_spin_rotation_vector; },
            "reference spin rotation vector (valid after a reference-particle push)"
        )
    ;
    register_push(py_SoftSolenoid);
    register_reverse(py_SoftSolenoid);

    py::class_<Source, elements::mixin::Named, elements::mixin::Thin> py_Source(me, "Source");
    py_Source
        .def("__repr__",
             [](Source const & src) {
                 return element_name(
                     src,
                     std::make_pair("distribution", src.m_distribution),
                     std::make_pair("openpmd_path", src.m_series_name),
                     std::make_pair("active_once", src.m_active_once),
                     std::make_pair("load_ref_particle", src.m_load_ref_particle),
                     std::make_pair("load_step", src.m_load_step),
                     std::make_pair("load_step_index", src.m_load_step_index)
                 );
             }
        )
        .def("to_dict",
            [](Source const & src) {
                // an unset option is None, so that the dict can be passed to the constructor
                ElementPropertyTypes load_step = py::none();
                if (src.m_load_step.has_value()) { load_step = src.m_load_step.value(); }
                ElementPropertyTypes load_step_index = py::none();
                if (src.m_load_step_index.has_value()) { load_step_index = src.m_load_step_index.value(); }

                return element_dict(
                    src,
                    std::make_pair("distribution", src.m_distribution),
                    std::make_pair("openpmd_path", src.m_series_name),
                    std::make_pair("active_once", src.m_active_once),
                    std::make_pair("load_ref_particle", src.m_load_ref_particle),
                    std::make_pair("load_step", load_step),
                    std::make_pair("load_step_index", load_step_index)
                );
            }
        )
        .def(py::init<
             std::string,
             std::string,
             bool,
             bool,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>
         >(),
             py::arg("distribution"),
             py::arg("openpmd_path"),
             py::arg("active_once") = Source::DEFAULT_active_once,
             py::arg("load_ref_particle") = Source::DEFAULT_load_ref_particle,
             py::arg("load_step") = Source::DEFAULT_load_step,
             py::arg("load_step_index") = Source::DEFAULT_load_step_index,
             py::arg("name") = py::none(),
             "A particle source."
        )
        .def_property("distribution",
            [](Source & src) { return src.m_distribution; },
            [](Source & src, std::string distribution) { src.m_distribution = distribution; },
            "Distribution type of particles in the source"
        )
        .def_property("series_name",
            [](Source & src) { return src.m_series_name; },
            [](Source & src, std::string series_name) { src.m_series_name = series_name; },
            "Path to openPMD series as accepted by openPMD_api.Series"
        )
        .def_property("active_once",
            [](Source & src) { return src.m_active_once; },
            [](Source & src, bool actice_once) { src.m_active_once = actice_once; },
            "Inject particles only for the first lattice period."
        )
        .def_property("load_ref_particle",
            [](Source & src) { return src.m_load_ref_particle; },
            [](Source & src, bool load_ref_particle) { src.m_load_ref_particle = load_ref_particle; },
            "Restore the reference particle from the species metadata of the openPMD file (particle tracking only)."
        )
        .def_property("load_step",
            [](Source & src) { return src.m_load_step; },
            [](Source & src, std::optional<int> load_step) { src.m_load_step = load_step; },
            "Which step (iteration) to load from the openPMD series: the ImpactX step at which "
            "the beam monitor wrote the beam, which is stored as the openPMD iteration in the "
            "file. Set at most one of load_step and load_step_index; if neither is set, the "
            "last step in the file is loaded."
        )
        .def_property("load_step_index",
            [](Source & src) { return src.m_load_step_index; },
            [](Source & src, std::optional<int> load_step_index) { src.m_load_step_index = load_step_index; },
            "Which step (iteration) to load from the openPMD series, by position in the file: "
            "0 is the first step and -1 the last, counting back from it as in Python "
            "(-2 is the second to last step). Set at most one of load_step and load_step_index; "
            "if neither is set, the last step in the file is loaded."
        )
    ;
    register_push(py_Source);
    register_reverse(py_Source);

    py::class_<Sol, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_Sol(me, "Sol");
    py_Sol
        .def("__repr__",
             [](Sol const & sol) {
                 return element_name(
                     sol,
                     std::make_pair("ks", sol.m_ks)
                 );
             }
        )
        .def("to_dict",
            [](Sol const & sol) {
                return element_dict(
                    sol,
                    std::make_pair("ks", sol.m_ks)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("ks"),
             py::arg("dx") = Sol::DEFAULT_dx,
             py::arg("dy") = Sol::DEFAULT_dy,
             py::arg("rotation") = Sol::DEFAULT_rotation_degree,
             py::arg("aperture_x") = Sol::DEFAULT_aperture_x,
             py::arg("aperture_y") = Sol::DEFAULT_aperture_y,
             py::arg("nslice") = Sol::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "An ideal hard-edge Solenoid magnet."
        )
        .def_property("ks",
            [](Sol & soft_sol) { return soft_sol.m_ks; },
            [](Sol & soft_sol, amrex::ParticleReal ks) { soft_sol.m_ks = ks; },
            "Solenoid strength in m^(-1) (MADX convention) in (magnetic field Bz in T) / (rigidity in T-m)"
        )
    ;
    register_push(py_Sol);
    register_reverse(py_Sol);

    py::class_<PRot, elements::mixin::Named, elements::mixin::Thin> py_PRot(me, "PRot");
    py_PRot
        .def("__repr__",
             [](PRot const & prot) {
                 return element_name(
                     prot,
                     std::make_pair("phi_in", prot.m_phi_in / elements::mixin::Alignment::degree2rad),
                     std::make_pair("phi_out", prot.m_phi_out / elements::mixin::Alignment::degree2rad)
                 );
             }
        )
        .def("to_dict",
            [](PRot const & prot, bool in_degrees) {
                if (in_degrees) {
                    return element_dict(
                        prot,
                        std::make_pair("phi_in", prot.m_phi_in / elements::mixin::Alignment::degree2rad),
                                                                     // once fixed, update src/python/impactx/extensions/KnownElementsList.py
                        std::make_pair("phi_out", prot.m_phi_out / elements::mixin::Alignment::degree2rad)
                    );
                } else {
                    // legacy: buggy radians instead of degrees
                    py::warnings::warn(
                        "Warning: PRot.to_dict() has a known bug, "
                        "returning phi_in and phi_out in radians than degrees. "
                        "Please use PRot.to_dict(in_degrees=True) instead.",
                        PyExc_RuntimeWarning,
                        2
                    );
                    return element_dict(
                        prot,
                        std::make_pair("phi_in", prot.m_phi_in),   // BUG: constructor is in degrees
                        std::make_pair("phi_out", prot.m_phi_out)  // BUG: constructor is in degrees
                    );
                }
            },
            py::arg("in_degrees") = false
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("phi_in"),
             py::arg("phi_out"),
             py::arg("name") = py::none(),
             "An exact pole-face rotation in the x-z plane. Both angles are in degrees."
        )
        .def_property("phi_in",
            [](PRot & prot) { return prot.m_phi_in; },
            [](PRot & prot, amrex::ParticleReal phi_in) { prot.m_phi_in = phi_in; },
            "angle of the reference particle with respect to the longitudinal (z) axis in the original frame in radian"
        )
        .def_property("phi_out",
            [](PRot & prot) { return prot.m_phi_out; },
            [](PRot & prot, amrex::ParticleReal phi_out) { prot.m_phi_out = phi_out; },
            "angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in radian"
        )
        /* BUG: this should be in degree
        .def_property("phi_in",
            [](PRot & prot) { return prot.m_phi_in / elements::mixin::Alignment::degree2rad; },
            [](PRot & prot, amrex::ParticleReal phi_in_deg) { prot.m_phi_in = phi_in_deg * elements::mixin::Alignment::degree2rad; },
            "angle of the reference particle with respect to the longitudinal (z) axis in the original frame in degrees"
        )
        .def_property("phi_out",
            [](PRot & prot) { return prot.m_phi_out / elements::mixin::Alignment::degree2rad; },
            [](PRot & prot, amrex::ParticleReal phi_out_deg) { prot.m_phi_out = phi_out_deg * elements::mixin::Alignment::degree2rad; },
            "angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in degrees"
        )
        */
    ;
    register_push(py_PRot);
    register_reverse(py_PRot);

    py::class_<SoftQuadrupole, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture> py_SoftQuadrupole(me, "SoftQuadrupole");
    py_SoftQuadrupole
        .def("__repr__",
             [](SoftQuadrupole const & soft_quad) {
                 return element_name(
                     soft_quad,
                     std::make_pair("gscale", soft_quad.m_gscale)
                 );
             }
        )
        .def("to_dict",
            [](SoftQuadrupole const & soft_quad) {
                return element_dict(
                    soft_quad,
                    std::make_pair("gscale", soft_quad.m_gscale),
                    std::make_pair("cos_coefficients", SoftQuadrupole::DynamicData::get(soft_quad.m_id)->cos.host_const()),
                    std::make_pair("sin_coefficients", SoftQuadrupole::DynamicData::get(soft_quad.m_id)->sin.host_const()),
                    std::make_pair("mapsteps", soft_quad.m_mapsteps)
                );
            }
        )
        .def(py::init<
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 std::vector<amrex::ParticleReal>,
                 std::vector<amrex::ParticleReal>,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 int,
                 int,
                 std::optional<std::string>
             >(),
             py::arg("ds"),
             py::arg("gscale"),
             py::arg("cos_coefficients"),
             py::arg("sin_coefficients"),
             py::arg("dx") = SoftQuadrupole::DEFAULT_dx,
             py::arg("dy") = SoftQuadrupole::DEFAULT_dy,
             py::arg("rotation") = SoftQuadrupole::DEFAULT_rotation_degree,
             py::arg("aperture_x") = SoftQuadrupole::DEFAULT_aperture_x,
             py::arg("aperture_y") = SoftQuadrupole::DEFAULT_aperture_y,
             py::arg("mapsteps") = SoftQuadrupole::DEFAULT_mapsteps,
             py::arg("nslice") = SoftQuadrupole::DEFAULT_nslice,
             py::arg("name") = py::none(),
             "A soft-edge quadrupole."
        )
        .def_property("gscale",
            [](SoftQuadrupole & soft_quad) { return soft_quad.m_gscale; },
            [](SoftQuadrupole & soft_quad, amrex::ParticleReal gscale) { soft_quad.m_gscale = gscale; },
            "Scaling factor for on-axis field gradient in inverse meters"
        )
        // TODO cos_coefficients
        // TODO sin_coefficients
        .def_property("mapsteps",
            [](SoftQuadrupole & soft_quad) { return soft_quad.m_mapsteps; },
            [](SoftQuadrupole & soft_quad, int mapsteps) {
                if (mapsteps <= 0)
                    throw std::invalid_argument("SoftQuadrupole: mapsteps must be > 0!");
                soft_quad.m_mapsteps = mapsteps;
            },
            "number of integration steps per slice used for map and reference particle push in applied fields"
        )
        .def_property_readonly("map",
            [](SoftQuadrupole const & soft_quad) { return soft_quad.m_map; },
            "linearized transport map around the reference particle (valid after a reference-particle push)"
        )
        .def_property_readonly("spin_coupling",
            [](SoftQuadrupole const & soft_quad) { return soft_quad.m_spin_coupling; },
            "linearized spin-orbit coupling matrix (valid after a reference-particle push)"
        )
    ;
    register_push(py_SoftQuadrupole);
    register_reverse(py_SoftQuadrupole);

    py::class_<ThinDipole, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_ThinDipole(me, "ThinDipole");
    py_ThinDipole
        .def("__repr__",
             [](ThinDipole const & thin_dp) {
                 return element_name(
                     thin_dp,
                     std::make_pair("theta", thin_dp.m_theta / ThinDipole::degree2rad),
                     std::make_pair("rc", thin_dp.m_rc)
                 );
             }
        )
        .def("to_dict",
            [](ThinDipole const & thin_dp, bool in_degrees) {
                if (in_degrees) {
                    return element_dict(
                        thin_dp,
                        std::make_pair("theta", thin_dp.m_theta / ThinDipole::degree2rad),
                                                                     // once fixed, update src/python/impactx/extensions/KnownElementsList.py
                        std::make_pair("rc", thin_dp.m_rc)
                    );
                } else {
                    // legacy: buggy radians instead of degrees
                    py::warnings::warn(
                        "Warning: ThinDipole.to_dict() has a known bug, "
                        "returning theta in radians than degrees. "
                        "Please use ThinDipole.to_dict(in_degrees=True) instead.",
                        PyExc_RuntimeWarning,
                        2
                    );
                    return element_dict(
                        thin_dp,
                        std::make_pair("theta", thin_dp.m_theta),  // BUG: constructor is in degrees
                        std::make_pair("rc", thin_dp.m_rc)
                    );
                }
            },
            py::arg("in_degrees") = false
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("theta"),
             py::arg("rc"),
             py::arg("dx") = ThinDipole::DEFAULT_dx,
             py::arg("dy") = ThinDipole::DEFAULT_dy,
             py::arg("rotation") = ThinDipole::DEFAULT_rotation_degree,
             py::arg("aperture_x") = ThinDipole::DEFAULT_aperture_x,
             py::arg("aperture_y") = ThinDipole::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             "A thin kick model of a dipole bend."
        )
        .def_property("theta",
            [](ThinDipole & thin_dp) { return thin_dp.m_theta; },
            [](ThinDipole & thin_dp, amrex::ParticleReal theta) { thin_dp.m_theta = theta; },
            "Bend angle (radian)"
        )
        /* BUG: this should be in degree
        .def_property("theta",
            [](ThinDipole & thin_dp) { return thin_dp.m_theta / ThinDipole::degree2rad; },
            [](ThinDipole & thin_dp, amrex::ParticleReal theta_deg) {
                thin_dp.m_theta = theta_deg * ThinDipole::degree2rad;
            },
            "Bend angle (degrees)"
        )
        .def_property("rc",
            [](ThinDipole & thin_dp) { return thin_dp.m_rc; },
            [](ThinDipole & thin_dp, amrex::ParticleReal rc) { thin_dp.m_rc = rc; },
            "Effective curvature radius (meters)"
        )
        */
    ;
    register_push(py_ThinDipole);
    register_reverse(py_ThinDipole);

    py::class_<TaperedPL, elements::mixin::Named, elements::mixin::Thin, elements::mixin::Alignment, elements::mixin::PipeAperture> py_TaperedPL(me, "TaperedPL");
    py_TaperedPL
        .def("__repr__",
             [](TaperedPL const & taperedpl) {
                 return element_name(
                     taperedpl,
                     std::make_pair("k", taperedpl.m_k),
                     std::make_pair("taper", taperedpl.m_taper)
                 );
             }
        )
        .def("to_dict",
            [](TaperedPL const & taperedpl) {
                return element_dict(
                    taperedpl,
                    std::make_pair("k", taperedpl.m_k),
                    std::make_pair("taper", taperedpl.m_taper),
                    std::make_pair("unit", taperedpl.m_unit)
                );
            }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("k"),
             py::arg("taper"),
             py::arg("unit") = TaperedPL::DEFAULT_unit,
             py::arg("dx") = TaperedPL::DEFAULT_dx,
             py::arg("dy") = TaperedPL::DEFAULT_dy,
             py::arg("rotation") = TaperedPL::DEFAULT_rotation_degree,
             py::arg("aperture_x") = TaperedPL::DEFAULT_aperture_x,
             py::arg("aperture_y") = TaperedPL::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             R"doc(A thin nonlinear plasma lens with transverse (horizontal) taper

             .. math::

                B_x = g \left( y + \frac{xy}{D_x} \right), \quad \quad B_y = -g \left(x + \frac{x^2 + y^2}{2 D_x} \right)

             where :math:`g` is the (linear) field gradient in T/m and :math:`D_x` is the targeted horizontal dispersion in m.
             )doc"
        )
        .def_property("k",
            [](TaperedPL & taperedpl) { return taperedpl.m_k; },
            [](TaperedPL & taperedpl, amrex::ParticleReal k) { taperedpl.m_k = k; },
            "integrated focusing strength in m^(-1) (if unit = 0) or integrated focusing strength in T (if unit = 1)"
        )
        .def_property("taper",
            [](TaperedPL & taperedpl) { return taperedpl.m_taper; },
            [](TaperedPL & taperedpl, amrex::ParticleReal taper) { taperedpl.m_taper = taper; },
            "horizontal taper parameter in m^(-1) = 1 / (target horizontal dispersion in m)"
        )
        .def_property("unit",
            [](TaperedPL & taperedpl) { return taperedpl.m_unit; },
            [](TaperedPL & taperedpl, int unit) { taperedpl.m_unit = unit; },
            "specification of units for plasma lens focusing strength"
        )
    ;
    register_push(py_TaperedPL);
    register_reverse(py_TaperedPL);

    py::class_<LinearMap, elements::mixin::Named, elements::mixin::Alignment, elements::mixin::PipeAperture> py_LinearMap(me, "LinearMap");
    py_LinearMap
        .def("__repr__",
             [](LinearMap const & linearmap) {
                 return element_name(linearmap);
             }
        )
        .def("to_dict",
            [](LinearMap const & linearmap) {
                return element_dict(
                    linearmap,
                    std::make_pair("ds", linearmap.m_ds),
                    std::make_pair("R", linearmap.m_transport_map)
                );
            }
        )
        .def(py::init<
                Map6x6,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("R"),
             py::arg("ds") = LinearMap::DEFAULT_ds,
             py::arg("dx") = LinearMap::DEFAULT_dx,
             py::arg("dy") = LinearMap::DEFAULT_dy,
             py::arg("rotation") = LinearMap::DEFAULT_rotation_degree,
             py::arg("aperture_x") = LinearMap::DEFAULT_aperture_x,
             py::arg("aperture_y") = LinearMap::DEFAULT_aperture_y,
             py::arg("name") = py::none(),
             "(A user-provided linear map, represented as a 6x6 transport matrix.)"
        )
        .def_property("R",
            [](LinearMap & linearmap) { return linearmap.m_transport_map; },
            [](LinearMap & linearmap, Map6x6 R) { linearmap.m_transport_map = R; },
            "linear map as a 6x6 transport matrix"
        )
        .def_property("ds",
            [](LinearMap & linearmap) { return linearmap.m_ds; },
            [](LinearMap & linearmap, amrex::ParticleReal ds) { linearmap.m_ds = ds; },
            "segment length in m"
        )
        .def_property_readonly("nslice",
            [](LinearMap & linearmap) { return linearmap.nslice(); },
            "one, because we do not support slicing of this element"
        )
        .def_property_readonly("symplectic", &LinearMap::symplectic,
            "Check if the transport map is symplectic.\n\n"
            "A matrix R is symplectic if R^T J R = J, where J is the\n"
            "standard 6x6 skew-symmetric symplectic form (also called Omega)."
        )
     ;
     register_push(py_LinearMap);
     register_reverse(py_LinearMap);

    py::class_<SpinMap, elements::mixin::Named, elements::mixin::Alignment> py_SpinMap(me, "SpinMap");
    py_SpinMap
        .def("__repr__",
             [](SpinMap const & spinmap) {
                 return element_name(spinmap);
             }
        )
        .def("to_dict",
            [](SpinMap const & spinmap) {
                return element_dict(
                    spinmap,
                    std::make_pair("ds", spinmap.m_ds),
                    std::make_pair("v", spinmap.m_spin_rotation_vector),
                    std::make_pair("A", spinmap.m_spin_orbit_coupling)
                );
            }
        )
        .def(py::init<
                Vector3,
                Map3x6,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::optional<std::string>
             >(),
             py::arg("v"),
             py::arg("A"),
             py::arg("ds") = SpinMap::DEFAULT_ds,
             py::arg("dx") = SpinMap::DEFAULT_dx,
             py::arg("dy") = SpinMap::DEFAULT_dy,
             py::arg("rotation") = SpinMap::DEFAULT_rotation_degree,
             py::arg("name") = py::none(),
             "(A user-provided spin map, represented as a 3-vector and a 3x6 coupling matrix.)"
        )
        .def_property("v",
            [](SpinMap & spinmap) { return spinmap.m_spin_rotation_vector; },
            [](SpinMap & spinmap, Vector3 v) { spinmap.m_spin_rotation_vector = v; },
            "design axis-angle generator of spin rotation as a 3x1 vector"
        )
        .def_property("A",
            [](SpinMap & spinmap) { return spinmap.m_spin_orbit_coupling; },
            [](SpinMap & spinmap, Map3x6 A) { spinmap.m_spin_orbit_coupling = A; },
            "spin-orbit coupling generator of rotation as a 3x6 matrix"
        )
        .def_property("ds",
            [](SpinMap & spinmap) { return spinmap.m_ds; },
            [](SpinMap & spinmap, amrex::ParticleReal ds) { spinmap.m_ds = ds; },
            "segment length in m"
        )
        .def_property_readonly("nslice",
            [](SpinMap & spinmap) { return spinmap.nslice(); },
            "one, because we do not support slicing of this element"
        )
     ;
     register_push(py_SpinMap);
     register_reverse(py_SpinMap);


    // freestanding push function
    m.def("push", py::overload_cast<ImpactXParticleContainer &, elements::KnownElements &, int, int>(&push),
        py::arg("pc"), py::arg("element"), py::arg("step")=0, py::arg("period")=0,
        "Push a whole particle beam (incl. reference particle) through an element"
    );
    m.def("push", py::overload_cast<RefPart &, elements::KnownElements &>(&push),
        py::arg("ref"), py::arg("element"),
        "Push the reference particle through an element"
    );

    m.def("reverse", [](elements::KnownElements & el) {
            elements::reverse(el);
        },
        py::arg("element"),
        "Reverse an element in-place so that pushing particles through\n"
        "it reverses the effect of the original element."
    );

    // init lattice next, under the elements sub-module
    init_lattice(me);
}
