/* Copyright 2022-2026 The ImpactX Community
 *
 * Authors: Ryan Sandberg, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ReferenceParticle.H>
#include <particles/transformation/CoordinateTransformation.H>

namespace py = pybind11;
using namespace impactx;

void init_transformation(py::module& m)
{
    m.def("coordinate_transformation",
        &transformation::CoordinateTransformation,
        py::arg("pc"),
        py::arg("direction"),
        "Transform coordinates from fixed s to fixed t or vice versa, using the reference energy pc.ref.pt at the time of the call."
    );
}
