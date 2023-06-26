/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/elements/All.H>
#include <AMReX.H>

#include <vector>

namespace py = pybind11;
using namespace impactx;


void init_elements(py::module& m)
{
    py::module_ me = m.def_submodule(
        "elements",
        "Accelerator lattice elements in ImpactX"
    );

    using KnownElementsList = std::list<KnownElements>;
    py::class_<KnownElementsList> kel(me, "KnownElementsList");
    kel
        .def(py::init<>())
        .def(py::init<KnownElements>())
        .def(py::init([](py::list l){
            auto v = new KnownElementsList;
            for (auto const & handle : l)
                v->push_back(handle.cast<KnownElements>());
            return v;
        }))

        .def("append", [](KnownElementsList &v, KnownElements el) { v.emplace_back(el); },
             "Add a single element to the list.")

        .def("extend",
             [](KnownElementsList &v, KnownElementsList l) {
                 for (auto const & el : l)
                    v.push_back(el);
                 return v;
             },
             "Add a list of elements to the list.")
        .def("extend",
             [](KnownElementsList &v, py::list l) {
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
                 amrex::ParticleReal const,
                 amrex::ParticleReal const
             >(),
             py::arg("ds"), py::arg("nslice") = 1,
             "Mixin class for lattice elements with finite length."
        )
        .def_property_readonly("nslice", &elements::Thick::nslice)
        .def_property_readonly("ds", &elements::Thick::ds)
    ;

    py::class_<elements::Thin>(me, "Thin")
        .def(py::init<>(),
             "Mixin class for lattice elements with zero length."
        )
        .def_property_readonly("nslice", &elements::Thin::nslice)
        .def_property_readonly("ds", &elements::Thin::ds)
    ;

    // diagnostics

    py::class_<diagnostics::BeamMonitor, elements::Thin>(me, "BeamMonitor")
        .def(py::init<std::string, std::string, std::string>(),
             py::arg("name"), py::arg("backend")="default", py::arg("encoding")="g",
             "This element writes the particle beam out to openPMD data."
        )
    ;

    // beam optics

    py::class_<ChrDrift, elements::Thick>(me, "ChrDrift")
        .def(py::init<
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("nslice") = 1,
             "A Drift with chromatic effects included."
        )
    ;

    py::class_<ChrQuad, elements::Thick>(me, "ChrQuad")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const,
        int const>(),
             py::arg("ds"), py::arg("k"), py::arg("units") = 0, py::arg("nslice") = 1,
             "A Quadrupole magnet with chromatic effects included."
        )
    ;

    py::class_<ChrAcc, elements::Thick>(me, "ChrAcc")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("ez"), py::arg("bz"), py::arg("nslice") = 1,
             "A region of Uniform Acceleration, with chromatic effects included."
        )
    ;

    py::class_<ConstF, elements::Thick>(me, "ConstF")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("kx"), py::arg("ky"), py::arg("kt"), py::arg("nslice") = 1,
             "A linear Constant Focusing element."
        )
    ;

    py::class_<DipEdge, elements::Thin>(me, "DipEdge")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("psi"), py::arg("rc"), py::arg("g"), py::arg("K2"),
             "Edge focusing associated with bend entry or exit."
        )
    ;

    py::class_<Drift, elements::Thick>(me, "Drift")
        .def(py::init<
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("nslice") = 1,
             "A drift."
        )
    ;

    py::class_<ExactDrift, elements::Thick>(me, "ExactDrift")
        .def(py::init<
                amrex::ParticleReal const,
                int const >(),
             py::arg("ds"), py::arg("nslice") = 1,
             "A Drift using the exact nonlinear map."
        )
    ;

    py::class_<ExactSbend, elements::Thick>(me, "ExactSbend")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("phi"), py::arg("B") = 0.0, py::arg("nslice") = 1,
             "An ideal sector bend using the exact nonlinear map.  When B = 0, the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B."
        )
    ;

    py::class_<Multipole, elements::Thin>(me, "Multipole")
        .def(py::init<
                int const,
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("multiple"), py::arg("K_normal"), py::arg("K_skew"),
             "A general thin multipole element."
        )
    ;

    py::class_<None, elements::Thin>(me, "None")
        .def(py::init<>(),
             "This element does nothing."
        )
    ;

    py::class_<NonlinearLens, elements::Thin>(me, "NonlinearLens")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("knll"), py::arg("cnll"),
             "Single short segment of the nonlinear magnetic insert element."
        )
    ;

    py::class_<Programmable>(me, "Programmable")
        .def(py::init<
                 amrex::ParticleReal,
                 int>(),
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
              ) { p.m_push = new_hook; },
              "hook for push of whole container (pc, step)"
        )
        .def_property("beam_particles",
              [](Programmable & p) { return p.m_beam_particles; },
              [](Programmable & p,
                 std::function<void(ImpactXParticleContainer::iterator *, RefPart &)> new_hook
              ) { p.m_beam_particles = new_hook; },
              "hook for beam particles (pti, RefPart)"
        )
        .def_property("ref_particle",
              [](Programmable & p) { return p.m_ref_particle; },
              [](Programmable & p,
                 std::function<void(RefPart &)> new_hook
              ) { p.m_ref_particle = new_hook; },
              "hook for reference particle (RefPart)"
        )
    ;

    py::class_<Quad, elements::Thick>(me, "Quad")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("k"), py::arg("nslice") = 1,
             "A Quadrupole magnet."
        )
    ;

    py::class_<RFCavity, elements::Thick>(me, "RFCavity")
        .def(py::init<
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                amrex::ParticleReal,
                std::vector<amrex::ParticleReal>,
                std::vector<amrex::ParticleReal>,
                int,
                int
             >(),
             py::arg("ds"), py::arg("escale"), py::arg("freq"),
             py::arg("phase"), py::arg("cos_coefficients"), py::arg("sin_coefficients"),
             py::arg("mapsteps") = 1, py::arg("nslice") = 1,
             "An RF cavity."
        )
    ;

    py::class_<Sbend, elements::Thick>(me, "Sbend")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("rc"), py::arg("nslice") = 1,
             "An ideal sector bend."
        )
    ;

    py::class_<CFbend, elements::Thick>(me, "CFbend")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
		amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("rc"), py::arg("k"), py::arg("nslice") = 1,
             "An ideal combined function bend (sector bend with quadrupole component)."
        )    
    ;        
        
    py::class_<ShortRF, elements::Thin>(me, "ShortRF")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("V"), py::arg("k"),
             "A short RF cavity element at zero crossing for bunching."
        )
    ;

    py::class_<SoftSolenoid, elements::Thick>(me, "SoftSolenoid")
        .def(py::init<
                 amrex::ParticleReal const,
                 amrex::ParticleReal const,
                 std::vector<amrex::ParticleReal>,
                 std::vector<amrex::ParticleReal>,
         int const,
                 int const
             >(),
             py::arg("ds"), py::arg("bscale"),
             py::arg("cos_coefficients"), py::arg("sin_coefficients"),
             py::arg("mapsteps") = 1, py::arg("nslice") = 1,
             "A soft-edge solenoid."
        )
    ;

    py::class_<Sol, elements::Thick>(me, "Sol")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const,
                int const>(),
             py::arg("ds"), py::arg("ks"), py::arg("nslice") = 1,
             "An ideal hard-edge Solenoid magnet."
        )
    ;

    py::class_<PRot, elements::Thin>(me, "PRot")
        .def(py::init<
                amrex::ParticleReal const,
                amrex::ParticleReal const>(),
             py::arg("phi_in"), py::arg("phi_out"),
             "An exact pole-face rotation in the x-z plane. Both angles are in degrees."
        )
    ;

    py::class_<SoftQuadrupole, elements::Thick>(me, "SoftQuadrupole")
        .def(py::init<
                 amrex::ParticleReal,
                 amrex::ParticleReal,
                 std::vector<amrex::ParticleReal>,
                 std::vector<amrex::ParticleReal>,
                 int,
                 int
             >(),
             py::arg("ds"), py::arg("gscale"),
             py::arg("cos_coefficients"), py::arg("sin_coefficients"),
             py::arg("mapsteps") = 1, py::arg("nslice") = 1,
             "A soft-edge quadrupole."
        )
    ;

}
