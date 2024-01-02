/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/Push.H>
#include <particles/elements/All.H>
#include <AMReX.H>

#include <type_traits>
#include <utility>
#include <vector>

namespace py = pybind11;
using namespace impactx;


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
    ;
    register_beamoptics_push(py_BeamMonitor);

    // beam optics

    py::class_<Aperture, elements::Thin, elements::Alignment> py_Aperture(me, "Aperture");
    py_Aperture
        .def(py::init([](
                 amrex::ParticleReal xmax,
                 amrex::ParticleReal ymax,
                 std::string shape,
                 amrex::ParticleReal dx,
                 amrex::ParticleReal dy,
                 amrex::ParticleReal rotation_degree
             )
             {
                 if (shape != "rectangular" && shape != "elliptical")
                     throw std::runtime_error("shape must be \"rectangular\" or \"elliptical\"");

                 Aperture::Shape s = shape == "rectangular" ?
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
    ;
    register_beamoptics_push(py_Aperture);

    py::class_<ChrDrift, elements::Thick, elements::Alignment> py_ChrDrift(me, "ChrDrift");
    py_ChrDrift
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
             py::arg("units") = 0,
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("nslice") = 1,
             "A Quadrupole magnet with chromatic effects included."
        )
    ;
    register_beamoptics_push(py_ChrQuad);

    py::class_<ChrAcc, elements::Thick, elements::Alignment> py_ChrAcc(me, "ChrAcc");
    py_ChrAcc
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
    ;
    register_beamoptics_push(py_ChrAcc);

    py::class_<ConstF, elements::Thick, elements::Alignment> py_ConstF(me, "ConstF");
    py_ConstF
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
    ;
    register_beamoptics_push(py_DipEdge);

    py::class_<Drift, elements::Thick, elements::Alignment> py_Drift(me, "Drift");
    py_Drift
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
    ;
    register_beamoptics_push(py_ExactSbend);

    py::class_<Kicker, elements::Thin, elements::Alignment> py_Kicker(me, "Kicker");
    py_Kicker
        .def(py::init([](
                amrex::ParticleReal xkick,
                amrex::ParticleReal ykick,
                std::string const & units,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             )
             {
                 if (units != "dimensionless" && units != "T-m")
                     throw std::runtime_error(R"(units must be "dimensionless" or "T-m")");

                 Kicker::UnitSystem const u = units == "dimensionless" ?
                                            Kicker::UnitSystem::dimensionless :
                                            Kicker::UnitSystem::Tm;
                 return new Kicker(xkick, ykick, u);
             }),
             py::arg("xkick"),
             py::arg("ykick"),
             py::arg("units") = "dimensionless",
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             R"(A thin transverse kicker element. Kicks are for units "dimensionless" or in "T-m".)"
        )
    ;
    register_beamoptics_push(py_Kicker);

    py::class_<Multipole, elements::Thin, elements::Alignment> py_Multipole(me, "Multipole");
    py_Multipole
        .def(py::init<
                int,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("multiple"),
             py::arg("K_normal"),
             py::arg("K_skew"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             "A general thin multipole element."
        )
    ;
    register_beamoptics_push(py_Multipole);

    py::class_<None, elements::Thin> py_None(me, "None");
    py_None
        .def(py::init<>(),
             "This element does nothing."
        )
    ;
    register_beamoptics_push(py_None);

    py::class_<NonlinearLens, elements::Thin, elements::Alignment> py_NonlinearLens(me, "NonlinearLens");
    py_NonlinearLens
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
    ;
    register_beamoptics_push(py_NonlinearLens);

    py::class_<Programmable>(me, "Programmable", py::dynamic_attr())
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
    ;
    register_beamoptics_push(py_Quad);

    py::class_<RFCavity, elements::Thick, elements::Alignment> py_RFCavity(me, "RFCavity");
    py_RFCavity
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
    ;
    register_beamoptics_push(py_RFCavity);

    py::class_<Sbend, elements::Thick, elements::Alignment> py_Sbend(me, "Sbend");
    py_Sbend
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
    ;
    register_beamoptics_push(py_Sbend);

    py::class_<CFbend, elements::Thick, elements::Alignment> py_CFbend(me, "CFbend");
    py_CFbend
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
    ;
    register_beamoptics_push(py_CFbend);

    py::class_<Buncher, elements::Thin, elements::Alignment> py_Buncher(me, "Buncher");
    py_Buncher
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
    ;
    register_beamoptics_push(py_Buncher);

    py::class_<ShortRF, elements::Thin, elements::Alignment> py_ShortRF(me, "ShortRF");
    py_ShortRF
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
    ;
    register_beamoptics_push(py_ShortRF);

    py::class_<SoftSolenoid, elements::Thick, elements::Alignment> py_SoftSolenoid(me, "SoftSolenoid");
    py_SoftSolenoid
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
             py::arg("bscale"),
             py::arg("cos_coefficients"),
             py::arg("sin_coefficients"),
             py::arg("dx") = 0,
             py::arg("dy") = 0,
             py::arg("rotation") = 0,
             py::arg("mapsteps") = 1,
             py::arg("nslice") = 1,
             "A soft-edge solenoid."
        )
    ;
    register_beamoptics_push(py_SoftSolenoid);

    py::class_<Sol, elements::Thick, elements::Alignment> py_Sol(me, "Sol");
    py_Sol
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
    ;
    register_beamoptics_push(py_Sol);

    py::class_<PRot, elements::Thin> py_PRot(me, "PRot");
    py_PRot
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal
             >(),
             py::arg("phi_in"),
             py::arg("phi_out"),
             "An exact pole-face rotation in the x-z plane. Both angles are in degrees."
        )
    ;
    register_beamoptics_push(py_PRot);

    py::class_<SoftQuadrupole, elements::Thick, elements::Alignment> py_SoftQuadrupole(me, "SoftQuadrupole");
    py_SoftQuadrupole
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
    ;
    register_beamoptics_push(py_SoftQuadrupole);

    py::class_<ThinDipole, elements::Thin, elements::Alignment> py_ThinDipole(me, "ThinDipole");
    py_ThinDipole
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
    ;
    register_beamoptics_push(py_ThinDipole);


    // free-standing push function
    m.def("push", &Push,
        py::arg("pc"), py::arg("element"), py::arg("step")=0,
        "Push particles through an element"
    );
}
