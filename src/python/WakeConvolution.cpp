#include "pyImpactX.H"

#include <particles/wakefields/WakeConvolution.H>
#include <particles/ImpactXParticleContainer.H> //Includes all necessary AMReX headers
#include <initialization/InitDistribution.H>

namespace py = pybind11;
using namespace impactx;

void init_wakeconvolution(py::module &m)
{
    py::module_ md = m.def_submodule("wakeconvolution");

    md.def("unit_step", &unit_step, "Step function");
    md.def("alpha", &alpha, "Alpha function");
    md.def("w_t_rf", &w_t_rf, "Transverse Resistive Wall Wake Function");
    md.def("w_l_rf", &w_l_rf, "Longitudinal Resistive Wall Wake Function");
    md.def("w_l_csr", &w_l_csr, "CSR Wake Function");
    //md.def("convolve_fft", &convolve_fft, "FFT Convolution");

    //Use a lambda to wrap the convolve_fft function
    md.def("convolve_fft", [](
        const std::vector<double>& beam_profile,
        const std::vector<double>& wake_func,
        double delta_t,
        std::vector<double>& result,
        int padding_factor)
    {
        convolve_fft(beam_profile, wake_func, delta_t, result, padding_factor);
    }, "FFT Convolution");
}
