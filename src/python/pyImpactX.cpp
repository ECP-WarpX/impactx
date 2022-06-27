/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <ImpactX.H>
#include <particles/elements/All.H>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace impactx;


// forward declarations of exposed classes
void init_distribution(py::module&);
void init_elements(py::module&);
void init_ImpactX(py::module&);
void init_impactxparticlecontainer(py::module&);
void init_refparticle(py::module&);

PYBIND11_MODULE(impactx_pybind, m) {
    // make sure AMReX types are known
    py::module::import("amrex");

    m.doc() = R"pbdoc(
            impactx_pybind
            --------------
            .. currentmodule:: impactx_pybind

            .. autosummary::
               :toctree: _generate
               ImpactX
               distribution
               elements
    )pbdoc";

    // note: order from parent to child classes
    init_distribution(m);
    init_elements(m);
    init_refparticle(m);
    init_impactxparticlecontainer(m);
    init_ImpactX(m);

    // API runtime version
    //   note PEP-440 syntax: x.y.zaN but x.y.z.devN
#ifdef PYIMPACTX_VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(PYIMPACTX_VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

    // authors
    m.attr("__author__") =
        "Axel Huebl, Chad Mitchell, Ryan Sandberg, Ji Qiang, et al.";

    // API runtime build-time feature variants
    // m.attr("variants") = impactx::getVariants();
    // TODO allow to query runtime versions of all dependencies

    // license SPDX identifier
    m.attr("__license__") = "BSD-3-Clause-LBNL";

    // TODO broken numpy if not at least v1.15.0: raise warning
    // auto numpy = py::module::import("numpy");
    // auto npversion = numpy.attr("__version__");
    // std::cout << "numpy version: " << py::str(npversion) << std::endl;
}
