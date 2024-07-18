/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/Push.H>
#include <particles/elements/All.H>
#include <AMReX.H>

#include <utility>
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
        T value;
        // TODO: if array do queryarr
        // bool const has_name = amrex::ParmParse(prefix).queryarr(name.c_str(), value);
        bool const has_name = amrex::ParmParse(prefix).query(name.c_str(), value);

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
    void register_beamoptics_push(T_PyClass & cl)
    {
        using Element = typename T_PyClass::type;  // py::class<T, options...>

        cl.def("push",
            [](Element & el, ImpactXParticleContainer & pc, int step) {
                el(pc, step);
            },
            py::arg("pc"), py::arg("step")=0,
            "Push first the reference particle, then all other particles."
        );
    }
}

void init_elements(py::module& m)
{
    py::module_ const me = m.def_submodule(
        "elements",
        "Accelerator lattice elements in ImpactX"
    );

    // mixin classes

    py::class_<elements::Thick>(me, "Thick")
        .def(py::init<
                 amrex::ParticleReal,
                 amrex::ParticleReal
             >(),
             py::arg("ds"),
             py::arg("nslice") = 1,
             "Mixin class for lattice elements with finite length."
        )
        .def_property("ds",
            [](elements::Thick & th) { return th.m_ds; },
            [](elements::Thick & th, amrex::ParticleReal ds) { th.m_ds = ds; },
            "segment length in m"
        )
        .def_property("nslice",
            [](elements::Thick & th) { return th.m_nslice; },
            [](elements::Thick & th, int nslice) { th.m_nslice = nslice; },
            "number of slices used for the application of space charge"
        )
    ;

    py::class_<elements::Thin>(me, "Thin")
        .def(py::init<>(),
             "Mixin class for lattice elements with zero length."
        )
        .def_property_readonly("ds",
            &elements::Thin::ds,
            "segment length in m"
        )
        .def_property_readonly("nslice",
            &elements::Thin::nslice,
            "number of slices used for the application of space charge"
        )
    ;

    py::class_<elements::Alignment>(me, "Alignment")
        .def(py::init<>(),
             "Mixin class for lattice elements with horizontal/vertical alignment errors."
        )
        .def_property("dx",
            [](elements::Alignment & a) { return a.dx(); },
            [](elements::Alignment & a, amrex::ParticleReal dx) { a.m_dx = dx; },
            "horizontal translation error in m"
        )
        .def_property("dy",
            [](elements::Alignment & a) { return a.dy(); },
            [](elements::Alignment & a, amrex::ParticleReal dy) { a.m_dy = dy; },
            "vertical translation error in m"
        )
        .def_property("rotation",
            [](elements::Alignment & a) { return a.rotation(); },
            [](elements::Alignment & a, amrex::ParticleReal rotation_degree)
            {
                a.m_rotation = rotation_degree * elements::Alignment::degree2rad;
            },
            "rotation error in the transverse plane in degree"
        )
    ;

    // diagnostics

    py::class_<diagnostics::BeamMonitor, elements::Thin> py_BeamMonitor(me, "BeamMonitor");
    py_BeamMonitor
        .def(py::init<std::string, std::string, std::string>(),
             py::arg("name"),
             py::arg("backend") = "default",
             py::arg("encoding") = "g",
             "This element writes the particle beam out to openPMD data."
        )
        .def_property_readonly("name",
            &diagnostics::BeamMonitor::series_name,
            "name of the series"
        )
        .def_property("nonlinear_lens_invariants",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<bool>(bm.series_name(), "nonlinear_lens_invariants"); },
            [](diagnostics::BeamMonitor & bm, bool nonlinear_lens_invariants) {
                amrex::ParmParse pp_element(bm.series_name());
                pp_element.add("nonlinear_lens_invariants", nonlinear_lens_invariants);
            },
            "Compute and output the invariants H and I within the nonlinear magnetic insert element"
        )
        .def_property("alpha",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.series_name(), "alpha"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real alpha) {
                amrex::ParmParse pp_element(bm.series_name());
                pp_element.add("alpha", alpha);
            },
            "Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.\n"
            "Horizontal and vertical values must be equal."
        )
        .def_property("beta",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.series_name(), "beta"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real beta) {
                amrex::ParmParse pp_element(bm.series_name());
                pp_element.add("beta", beta);
            },
            "Twiss beta (in meters) of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.\n"
            "Horizontal and vertical values must be equal."
        )
        .def_property("tn",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.series_name(), "tn"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real tn) {
                amrex::ParmParse pp_element(bm.series_name());
                pp_element.add("tn", tn);
            },
            "Dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I."
        )
        .def_property("cn",
            [](diagnostics::BeamMonitor & bm) { return detail::get_or_throw<amrex::Real>(bm.series_name(), "cn"); },
            [](diagnostics::BeamMonitor & bm, amrex::Real cn) {
                amrex::ParmParse pp_element(bm.series_name());
                pp_element.add("cn", cn);
            },
            "Scale factor (in meters^(1/2)) of the IOTA nonlinear magnetic insert element used for computing H and I."
        )
    ;

    register_beamoptics_push(py_BeamMonitor);

    // beam optics

    py::class_<Aperture, elements::Thin, elements::Alignment> py_Aperture(me, "Aperture");
    py_Aperture
        .def("__repr__",
             [](Aperture const & /* ap */) {
                 return std::string("<impactx.elements.Aperture>");
             }
        )
        .def(py::init([](
                 amrex::ParticleReal xmax,
                 amrex::ParticleReal ymax,
                 std::string const & shape,
                 amrex::ParticleReal dx,
                 amrex::ParticleReal dy,
                 amrex::ParticleReal rotation_degree
             )
             {
                 if (shape != "rectangular" && shape != "elliptical")
                     throw std::runtime_error(R"(shape must be "rectangular" or "elliptical")");

                 Aperture::Shape const s = shape == "rectangular" ?
                     Aperture::Shape::rectangular :
                     Aperture::Shape::elliptical;
                 return new Aperture(xmax, ymax, s, dx, dy, rotation_degree);
             }),
             py::arg("xmax"),
             py::arg("ymax"),
             py::arg("shape") = "rectangular",
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             "A short collimator element applying a transverse aperture boundary."
        )
        .def_property("shape",
            [](Aperture & ap)
            {
                switch (ap.m_shape)
                {
                    case Aperture::Shape::rectangular :  // default
                        return "rectangular";
                    case Aperture::Shape::elliptical :
                        return "elliptical";
                    default:
                        throw std::runtime_error("Unknown shape");
                }
            },
            [](Aperture & ap, std::string const & shape)
            {
                if (shape != "rectangular" && shape != "elliptical")
                    throw std::runtime_error(R"(shape must be "rectangular" or "elliptical")");

                ap.m_shape = shape == "rectangular" ?
                    Aperture::Shape::rectangular :
                    Aperture::Shape::elliptical;
            },
            "aperture type (rectangular, elliptical)"
        )
        .def_property("xmax",
            [](Aperture & ap) { return ap.m_xmax; },
            [](Aperture & ap, amrex::ParticleReal xmax) { ap.m_xmax = xmax; },
            "maximum horizontal coordinate"
        )
        .def_property("ymax",
            [](Aperture & ap) { return ap.m_ymax; },
            [](Aperture & ap, amrex::ParticleReal ymax) { ap.m_ymax = ymax; },
            "maximum vertical coordinate"
        )
    ;
    register_beamoptics_push(py_Aperture);

    py::class_<ChrDrift, elements::Thick, elements::Alignment> py_ChrDrift(me, "ChrDrift");
    py_ChrDrift
        .def("__repr__",
             [](ChrDrift const & chr_drift) {
                 std::string r = "<impactx.elements.ChrDrift (ds=";
                 r.append(std::to_string(chr_drift.ds()))
                         .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "A Drift with chromatic effects included."
        )
    ;
    register_beamoptics_push(py_ChrDrift);

    py::class_<ChrQuad, elements::Thick, elements::Alignment> py_ChrQuad(me, "ChrQuad");
    py_ChrQuad
        .def("__repr__",
             [](ChrQuad const & chr_quad) {
                 std::string r = "<impactx.elements.ChrQuad (ds=";
                 r.append(std::to_string(chr_quad.ds()))
                  .append(", k=")
                  .append(std::to_string(chr_quad.m_k))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("k"),
             py::arg("unit") = 0,
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
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
    register_beamoptics_push(py_ChrQuad);

    py::class_<ChrPlasmaLens, elements::Thick, elements::Alignment> py_ChrPlasmaLens(me, "ChrPlasmaLens");
    py_ChrPlasmaLens
        .def("__repr__",
             [](ChrPlasmaLens const & chr_pl_lens) {
                 std::string r = "<impactx.elements.ChrPlasmaLens (ds=";
                 r.append(std::to_string(chr_pl_lens.ds()))
                  .append(", k=")
                  .append(std::to_string(chr_pl_lens.m_k))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("k"),
             py::arg("unit") = 0,
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "An active Plasma Lens with chromatic effects included."
        )
        .def_property("k",
            [](ChrQuad & cq) { return cq.m_k; },
            [](ChrQuad & cq, amrex::ParticleReal k) { cq.m_k = k; },
            "focusing strength in 1/m^2 (or T/m)"
        )
        .def_property("unit",
            [](ChrQuad & cq) { return cq.m_unit; },
            [](ChrQuad & cq, int unit) { cq.m_unit = unit; },
            "unit specification for focusing strength"
        )
    ;
    register_beamoptics_push(py_ChrPlasmaLens);

    py::class_<ChrAcc, elements::Thick, elements::Alignment> py_ChrAcc(me, "ChrAcc");
    py_ChrAcc
        .def("__repr__",
             [](ChrAcc const & chr_acc) {
                 std::string r = "<impactx.elements.ChrAcc (ds=";
                 r.append(std::to_string(chr_acc.ds()))
                  .append(", ez=")
                  .append(std::to_string(chr_acc.m_ez))
                  .append(", bz=")
                  .append(std::to_string(chr_acc.m_bz))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("ez"),
             py::arg("bz"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
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
    register_beamoptics_push(py_ChrAcc);

    py::class_<ConstF, elements::Thick, elements::Alignment> py_ConstF(me, "ConstF");
    py_ConstF
        .def("__repr__",
             [](ConstF const & constf) {
                 std::string r = "<impactx.elements.ConstF (ds=";
                 r.append(std::to_string(constf.ds()))
                  .append(", kx=")
                  .append(std::to_string(constf.m_kx))
                  .append(", ky=")
                  .append(std::to_string(constf.m_ky))
                  .append(", kt=")
                  .append(std::to_string(constf.m_kt))
                  .append(")>");
                 return r;
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
                int
             >(),
             py::arg("ds"),
             py::arg("kx"),
             py::arg("ky"),
             py::arg("kt"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
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
    register_beamoptics_push(py_ConstF);

    py::class_<DipEdge, elements::Thin, elements::Alignment> py_DipEdge(me, "DipEdge");
    py_DipEdge
        .def("__repr__",
             [](DipEdge const & dip_edge) {
                 std::string r = "<impactx.elements.DipEdge (psi=";
                 r.append(std::to_string(dip_edge.m_psi))
                  .append(", rc=")
                  .append(std::to_string(dip_edge.m_rc))
                  .append(", g=")
                  .append(std::to_string(dip_edge.m_g))
                  .append(", K2=")
                  .append(std::to_string(dip_edge.m_K2))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("psi"),
             py::arg("rc"),
             py::arg("g"),
             py::arg("K2"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
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
        .def_property("K2",
            [](DipEdge & dip_edge) { return dip_edge.m_K2; },
            [](DipEdge & dip_edge, amrex::ParticleReal K2) { dip_edge.m_K2 = K2; },
            "Fringe field integral (unitless)"
        )
    ;
    register_beamoptics_push(py_DipEdge);

    py::class_<Drift, elements::Thick, elements::Alignment> py_Drift(me, "Drift");
    py_Drift
        .def("__repr__",
             [](Drift const & drift) {
                 std::string r = "<impactx.elements.Drift (ds=";
                 r.append(std::to_string(drift.ds()))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "A drift."
        )
    ;
    register_beamoptics_push(py_Drift);

    py::class_<ExactDrift, elements::Thick, elements::Alignment> py_ExactDrift(me, "ExactDrift");
    py_ExactDrift
        .def("__repr__",
             [](ExactDrift const & exact_drift) {
                 std::string r = "<impactx.elements.ExactDrift (ds=";
                 r.append(std::to_string(exact_drift.ds()))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "A Drift using the exact nonlinear map."
        )
    ;
    register_beamoptics_push(py_ExactDrift);

    py::class_<ExactSbend, elements::Thick, elements::Alignment> py_ExactSbend(me, "ExactSbend");
    py_ExactSbend
        .def("__repr__",
             [](ExactSbend const & exact_sbend) {
                 std::string r = "<impactx.elements.ExactSbend (ds=";
                 r.append(std::to_string(exact_sbend.ds()))
                  .append(", phi=")
                  .append(std::to_string(exact_sbend.m_phi))
                  .append(", B=")
                  .append(std::to_string(exact_sbend.m_B))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("phi"),
             py::arg("B") = 0.0,
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "An ideal sector bend using the exact nonlinear map.  When B = 0, the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B."
        )
        .def_property("phi",
            [](ExactSbend & exact_sbend) { return exact_sbend.m_phi; },
            [](ExactSbend & exact_sbend, amrex::ParticleReal phi) { exact_sbend.m_phi = phi; },
            "Bend angle in degrees"
        )
        .def_property("B",
            [](ExactSbend & exact_sbend) { return exact_sbend.m_B; },
            [](ExactSbend & exact_sbend, amrex::ParticleReal B) { exact_sbend.m_B = B; },
            "Magnetic field in Tesla; when B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B"
        )
    ;
    register_beamoptics_push(py_ExactSbend);

    py::class_<Kicker, elements::Thin, elements::Alignment> py_Kicker(me, "Kicker");
    py_Kicker
        .def("__repr__",
             [](Kicker const & kicker) {
                 std::string r = "<impactx.elements.Kicker (xkick=";
                 r.append(std::to_string(kicker.m_xkick))
                  .append(", ykick=")
                  .append(std::to_string(kicker.m_ykick))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init([](
                amrex::ParticleReal xkick,
                amrex::ParticleReal ykick,
                std::string const & unit,
                amrex::ParticleReal dx,
                amrex::ParticleReal dy,
                amrex::ParticleReal rotation_degree
             )
             {
                 if (unit != "dimensionless" && unit != "T-m")
                     throw std::runtime_error(R"(unit must be "dimensionless" or "T-m")");

                 Kicker::UnitSystem const u = unit == "dimensionless" ?
                                            Kicker::UnitSystem::dimensionless :
                                            Kicker::UnitSystem::Tm;
                 return new Kicker(xkick, ykick, u, dx, dy, rotation_degree);
             }),
             py::arg("xkick"),
             py::arg("ykick"),
             py::arg("unit") = "dimensionless",
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
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
    register_beamoptics_push(py_Kicker);

    py::class_<Multipole, elements::Thin, elements::Alignment> py_Multipole(me, "Multipole");
    py_Multipole
        .def("__repr__",
             [](Multipole const & multipole) {
                 std::string r = "<impactx.elements.Multipole (multipole=";
                 r.append(std::to_string(multipole.m_multipole))
                  .append(", K_normal=")
                  .append(std::to_string(multipole.m_Kn))
                  .append(", K_skew=")
                  .append(std::to_string(multipole.m_Ks))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("multipole"),
             py::arg("K_normal"),
             py::arg("K_skew"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             "A general thin multipole element."
        )
        .def_property("multipole",
            [](Multipole & multipole) { return multipole.m_multipole; },
            [](Multipole & multipole, amrex::ParticleReal multipole_index) { multipole.m_multipole = multipole_index; },
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
    register_beamoptics_push(py_Multipole);

    py::class_<Empty, elements::Thin> py_None(me, "Empty");
    py_None
        .def("__repr__",
             [](Empty const & /* none */) {
                 return std::string("<impactx.elements.Empty>");
             }
        )
        .def(py::init<>(),
             "This element does nothing."
        )
    ;
    register_beamoptics_push(py_None);

    py::class_<NonlinearLens, elements::Thin, elements::Alignment> py_NonlinearLens(me, "NonlinearLens");
    py_NonlinearLens
        .def("__repr__",
             [](NonlinearLens const & nl) {
                 std::string r = "<impactx.elements.NonlinearLens (knll=";
                 r.append(std::to_string(nl.m_knll))
                  .append(", cnll=")
                  .append(std::to_string(nl.m_cnll))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("knll"),
             py::arg("cnll"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
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
    register_beamoptics_push(py_NonlinearLens);

    py::class_<Programmable>(me, "Programmable", py::dynamic_attr())
        .def("__repr__",
             [](Programmable const & prg) {
                 std::string r = "<impactx.elements.Programmable (ds=";
                 r.append(std::to_string(prg.ds()))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                 amrex::ParticleReal,
                 int
             >(),
             py::arg("ds") = 0.0, py::arg("nslice") = 1,
             "A programmable beam optics element."
        )
        .def_property("nslice",
            [](Programmable & p) { return p.nslice(); },
            [](Programmable & p, int nslice) { p.m_nslice = nslice; }
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
                 std::function<void(ImpactXParticleContainer *, int)> new_hook
              ) { p.m_push = std::move(new_hook); },
              "hook for push of whole container (pc, step)"
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

    py::class_<Quad, elements::Thick, elements::Alignment> py_Quad(me, "Quad");
    py_Quad
        .def("__repr__",
             [](Quad const & quad) {
                 std::string r = "<impactx.elements.Quad (ds=";
                 r.append(std::to_string(quad.ds()))
                  .append(", k=")
                  .append(std::to_string(quad.m_k))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("k"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "A Quadrupole magnet."
        )
        .def_property("k",
            [](Quad & quad) { return quad.m_k; },
            [](Quad & quad, amrex::ParticleReal k) { quad.m_k = k; },
            "Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing"
        )
    ;
    register_beamoptics_push(py_Quad);

    py::class_<RFCavity, elements::Thick, elements::Alignment> py_RFCavity(me, "RFCavity");
    py_RFCavity
        .def("__repr__",
             [](RFCavity const & rfc) {
                 std::string r = "<impactx.elements.RFCavity (ds=";
                 r.append(std::to_string(rfc.ds()))
                  .append(", escale=")
                  .append(std::to_string(rfc.m_escale))
                  .append(", freq=")
                  .append(std::to_string(rfc.m_freq))
                  .append(", phase=")
                  .append(std::to_string(rfc.m_phase))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::vector<amrex::ParticleReal>,
                std::vector<amrex::ParticleReal>,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                int
             >(),
             py::arg("ds"),
             py::arg("escale"),
             py::arg("freq"),
             py::arg("phase"),
             py::arg("cos_coefficients"),
             py::arg("sin_coefficients"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("mapsteps") = 1,
             py::arg("nslice") = 1,
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
            [](RFCavity & rfc, int mapsteps) { rfc.m_mapsteps = mapsteps; },
            "number of integration steps per slice used for map and reference particle push in applied fields"
        )
    ;
    register_beamoptics_push(py_RFCavity);

    py::class_<Sbend, elements::Thick, elements::Alignment> py_Sbend(me, "Sbend");
    py_Sbend
        .def("__repr__",
             [](Sbend const & sbend) {
                 std::string r = "<impactx.elements.Sbend (ds=";
                 r.append(std::to_string(sbend.ds()))
                  .append(", rc=")
                  .append(std::to_string(sbend.m_rc))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("rc"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "An ideal sector bend."
        )
        .def_property("rc",
            [](Sbend & sbend) { return sbend.m_rc; },
            [](Sbend & sbend, amrex::ParticleReal rc) { sbend.m_rc = rc; },
            "Radius of curvature in m"
        )
    ;
    register_beamoptics_push(py_Sbend);

    py::class_<CFbend, elements::Thick, elements::Alignment> py_CFbend(me, "CFbend");
    py_CFbend
        .def("__repr__",
             [](CFbend const & cfbend) {
                 std::string r = "<impactx.elements.CFbend (ds=";
                 r.append(std::to_string(cfbend.ds()))
                  .append(", rc=")
                  .append(std::to_string(cfbend.m_rc))
                  .append(", k=")
                  .append(std::to_string(cfbend.m_k))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("rc"),
             py::arg("k"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
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
    register_beamoptics_push(py_CFbend);

    py::class_<Buncher, elements::Thin, elements::Alignment> py_Buncher(me, "Buncher");
    py_Buncher
        .def("__repr__",
             [](Buncher const & buncher) {
                 std::string r = "<impactx.elements.Buncher (V=";
                 r.append(std::to_string(buncher.m_V))
                  .append(", k=")
                  .append(std::to_string(buncher.m_k))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("V"),
             py::arg("k"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
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
    register_beamoptics_push(py_Buncher);

    py::class_<ShortRF, elements::Thin, elements::Alignment> py_ShortRF(me, "ShortRF");
    py_ShortRF
        .def("__repr__",
             [](ShortRF const & short_rf) {
                 std::string r = "<impactx.elements.ShortRF (V=";
                 r.append(std::to_string(short_rf.m_V))
                  .append(", freq=")
                  .append(std::to_string(short_rf.m_freq))
                  .append(", phase=")
                  .append(std::to_string(short_rf.m_phase))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("V"),
             py::arg("freq"),
             py::arg("phase") = -90.0,
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
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
    register_beamoptics_push(py_ShortRF);

    py::class_<SoftSolenoid, elements::Thick, elements::Alignment> py_SoftSolenoid(me, "SoftSolenoid");
    py_SoftSolenoid
        .def("__repr__",
             [](SoftSolenoid const & soft_sol) {
                 std::string r = "<impactx.elements.SoftSolenoid (ds=";
                 r.append(std::to_string(soft_sol.ds()))
                  .append(", bscale=")
                  .append(std::to_string(soft_sol.m_bscale))
                  .append(")>");
                 return r;
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
                 int,
                 int
             >(),
             py::arg("ds"),
             py::arg("bscale"),
             py::arg("cos_coefficients"),
             py::arg("sin_coefficients"),
             py::arg("unit") = 0,
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("mapsteps") = 1,
             py::arg("nslice") = 1,
             "A soft-edge solenoid."
        )
        .def_property("bscale",
            [](SoftSolenoid & soft_sol) { return soft_sol.m_bscale; },
            [](SoftSolenoid & soft_sol, amrex::ParticleReal bscale) { soft_sol.m_bscale = bscale; },
            "Scaling factor for on-axis magnetic field Bz in inverse meters (if unit = 0) or magnetic field Bz in T (SI units, if unit = 1)"
        )
        // TODO cos_coefficients
        // TODO sin_coefficients
        .def_property("unit",
            [](SoftSolenoid & soft_sol) { return soft_sol.m_unit; },
            [](SoftSolenoid & soft_sol, amrex::ParticleReal unit) { soft_sol.m_unit = unit; },
            "specification of units for scaling of the on-axis longitudinal magnetic field"
        )
        .def_property("mapsteps",
            [](SoftSolenoid & soft_sol) { return soft_sol.m_mapsteps; },
            [](SoftSolenoid & soft_sol, int mapsteps) { soft_sol.m_mapsteps = mapsteps; },
            "number of integration steps per slice used for map and reference particle push in applied fields"
        )
    ;
    register_beamoptics_push(py_SoftSolenoid);

    py::class_<Sol, elements::Thick, elements::Alignment> py_Sol(me, "Sol");
    py_Sol
        .def("__repr__",
             [](Sol const & soft_sol) {
                 std::string r = "<impactx.elements.Sol (ds=";
                 r.append(std::to_string(soft_sol.ds()))
                  .append(", ks=")
                  .append(std::to_string(soft_sol.m_ks))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                int
             >(),
             py::arg("ds"),
             py::arg("ks"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "An ideal hard-edge Solenoid magnet."
        )
        .def_property("ks",
            [](Sol & soft_sol) { return soft_sol.m_ks; },
            [](Sol & soft_sol, amrex::ParticleReal ks) { soft_sol.m_ks = ks; },
            "Solenoid strength in m^(-1) (MADX convention) in (magnetic field Bz in T) / (rigidity in T-m)"
        )
    ;
    register_beamoptics_push(py_Sol);

    py::class_<PRot, elements::Thin> py_PRot(me, "PRot");
    py_PRot
        .def("__repr__",
             [](PRot const & prot) {
                 std::string r = "<impactx.elements.PRot (phi_in=";
                 r.append(std::to_string(prot.m_phi_in))
                  .append(", phi_out=")
                  .append(std::to_string(prot.m_phi_out))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("phi_in"),
             py::arg("phi_out"),
             "An exact pole-face rotation in the x-z plane. Both angles are in degrees."
        )
        .def_property("phi_in",
            [](PRot & prot) { return prot.m_phi_in; },
            [](PRot & prot, amrex::ParticleReal phi_in) { prot.m_phi_in = phi_in; },
            "angle of the reference particle with respect to the longitudinal (z) axis in the original frame in degrees"
        )
        .def_property("phi_out",
            [](PRot & prot) { return prot.m_phi_out; },
            [](PRot & prot, amrex::ParticleReal phi_out) { prot.m_phi_out = phi_out; },
            "angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in degrees"
        )
    ;
    register_beamoptics_push(py_PRot);

    py::class_<SoftQuadrupole, elements::Thick, elements::Alignment> py_SoftQuadrupole(me, "SoftQuadrupole");
    py_SoftQuadrupole
        .def("__repr__",
             [](SoftQuadrupole const & soft_quad) {
                 std::string r = "<impactx.elements.SoftQuadrupole (ds=";
                 r.append(std::to_string(soft_quad.ds()))
                  .append(", gscale=")
                  .append(std::to_string(soft_quad.m_gscale))
                  .append(")>");
                 return r;
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
                 int,
                 int
             >(),
             py::arg("ds"),
             py::arg("gscale"),
             py::arg("cos_coefficients"),
             py::arg("sin_coefficients"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("mapsteps") = 1,
             py::arg("nslice") = 1,
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
            [](SoftQuadrupole & soft_quad, int mapsteps) { soft_quad.m_mapsteps = mapsteps; },
            "number of integration steps per slice used for map and reference particle push in applied fields"
        )
    ;
    register_beamoptics_push(py_SoftQuadrupole);

    py::class_<ThinDipole, elements::Thin, elements::Alignment> py_ThinDipole(me, "ThinDipole");
    py_ThinDipole
        .def("__repr__",
             [](ThinDipole const & thin_dp) {
                 std::string r = "<impactx.elements.ThinDipole (theta=";
                 r.append(std::to_string(thin_dp.m_theta))
                  .append(", rc=")
                  .append(std::to_string(thin_dp.m_rc))
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("theta"),
             py::arg("rc"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             "A thin kick model of a dipole bend."
        )
        .def_property("theta",
            [](ThinDipole & thin_dp) { return thin_dp.m_theta; },
            [](ThinDipole & thin_dp, amrex::ParticleReal theta) { thin_dp.m_theta = theta; },
            "Bend angle (degrees)"
        )
        .def_property("rc",
            [](ThinDipole & thin_dp) { return thin_dp.m_rc; },
            [](ThinDipole & thin_dp, amrex::ParticleReal rc) { thin_dp.m_rc = rc; },
            "Effective curvature radius (meters)"
        )
    ;
    register_beamoptics_push(py_ThinDipole);

    py::class_<TaperedPL, elements::Thin, elements::Alignment> py_TaperedPL(me, "TaperedPL");
    py_TaperedPL
        .def("__repr__",
             [](TaperedPL const & taperedpl) {
                 std::string r = "<impactx.elements.TaperedPL (taperedpl=";
                 r.append(std::to_string(taperedpl.m_k))
                  .append(", k=")
                  .append(std::to_string(taperedpl.m_taper))
                  .append(", taper=")
                  .append(")>");
                 return r;
             }
        )
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("k"),
             py::arg("taper"),
             py::arg("unit") = 0,
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
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
    register_beamoptics_push(py_TaperedPL);


    // free-standing push function
    m.def("push", &Push,
        py::arg("pc"), py::arg("element"), py::arg("step")=0,
        "Push particles through an element"
    );


    // all-element type list
    using KnownElementsList = std::list<KnownElements>;
    py::class_<KnownElementsList> kel(me, "KnownElementsList");
    kel
        .def(py::init<>())
        .def(py::init<KnownElements>())
        .def(py::init([](py::list const & l){
            auto v = new KnownElementsList;
            for (auto const & handle : l)
                v->push_back(handle.cast<KnownElements>());
            return v;
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

        .def("clear", &KnownElementsList::clear,
             "Clear the list to become empty.")
        .def("pop_back", &KnownElementsList::pop_back,
             "Return and remove the last element of the list.")
        .def("__len__", [](const KnownElementsList &v) { return v.size(); },
             "The length of the list.")
        .def("__iter__", [](KnownElementsList &v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>()) /* Keep list alive while iterator is used */
    ;
}
