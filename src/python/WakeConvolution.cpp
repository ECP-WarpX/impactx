/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/wakefields/WakeConvolution.H>
#include <particles/wakefields/ChargeBinning.H>
#include <particles/ImpactXParticleContainer.H>
#include <initialization/InitDistribution.H>

#include <pybind11/numpy.h>  //Include pybind11 numpy header

namespace py = pybind11;
using namespace impactx;

void init_wakeconvolution(py::module &m)
{
    py::module_ md = m.def_submodule("wakeconvolution");

    //Use a lambda to wrap the binning functions
    md.def("deposit_charge",
        &impactx::particles::wakefields::DepositCharge1D,
        "Deposit Charge Distribution Function"
    );

    md.def("derivative_charge",
        &impactx::particles::wakefields::DerivativeCharge1D,
        "Derivative of Charge Profile Function"
    );

    md.def("unit_step", &impactx::particles::wakefields::unit_step, "Step Function");
    md.def("alpha", &impactx::particles::wakefields::alpha, "Alpha Function");
    md.def("w_t_rf", &impactx::particles::wakefields::w_t_rf, "Transverse Resistive Wall Wake Function");
    md.def("w_l_rf", &impactx::particles::wakefields::w_l_rf, "Longitudinal Resistive Wall Wake Function");
    md.def("w_l_csr", &impactx::particles::wakefields::w_l_csr, "CSR Wake Function");

    // Use a lambda to wrap the convolution function
    md.def("convolve_fft",
        &impactx::particles::wakefields::convolve_fft,
        "FFT Convolution"
    );
}
