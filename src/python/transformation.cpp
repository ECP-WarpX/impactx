/* Copyright 2021-2023 The ImpactX Community
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
        "Transform coordinates from fixed s to fixed to or vice versa."
    );

    py::enum_<transformation::Direction>(m, "TransformationDirection")
        .value("to_fixed_s", transformation::Direction::to_fixed_s)
        .value("to_fixed_t", transformation::Direction::to_fixed_t)
        .export_values();
}
