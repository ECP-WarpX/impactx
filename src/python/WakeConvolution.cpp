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
    md.def("deposit_charge", [](ImpactXParticleContainer& pc, py::array_t<double> charge_distribution, int num_bins, double bin_min, double bin_size, bool is_unity_particle_weight) {
        double* dptr_data = charge_distribution.mutable_data();
        DepositCharge1D(pc, dptr_data, num_bins, bin_min, bin_size, is_unity_particle_weight);
    }, "Deposit Charge Distribution Function");

    md.def("derivative_charge", [](py::array_t<double> charge_distribution, py::array_t<double> slopes, int num_bins, double bin_size, bool GetNumberDensity) {
        double* charge_dist = charge_distribution.mutable_data();
        double* derivatives = slopes.mutable_data();
        DerivativeCharge1D(charge_dist, derivatives, num_bins, bin_size, GetNumberDensity);
    }, "Derivative of Charge Profile Function");

    md.def("unit_step", &unit_step, "Step Function");
    md.def("alpha", &alpha, "Alpha Function");
    md.def("w_t_rf", &w_t_rf, "Transverse Resistive Wall Wake Function");
    md.def("w_l_rf", &w_l_rf, "Longitudinal Resistive Wall Wake Function");
    md.def("w_l_csr", &w_l_csr, "CSR Wake Function");

    //Use a lambda to wrap the convolution function
    md.def("convolve_fft", [](
        py::array_t<double> beam_profile,
        py::array_t<double> wake_func,
        double delta_t,
        py::array_t<double> result,
        int padding_factor)
    {
        double* beam_p = beam_profile.mutable_data();
        double* wake_f = wake_func.mutable_data();
        double* res = result.mutable_data();
        int beam_profile_size = beam_profile.size();
        int wake_func_size = wake_func.size();
        return convolve_fft(beam_p, wake_f, beam_profile_size, wake_func_size, delta_t, res, padding_factor);
    }, "FFT Convolution");
}