#ifndef WAKE_CONVOLUTION_H
#define WAKE_CONVOLUTION_H

#include "particles/ImpactXParticleContainer.H"
#include "initialization/InitDistribution.H"


//Constants
constexpr double Z0 = 377;  //Free space impedance [Ohm]
constexpr double alpha_1 = 0.4648;

//Step Function
double unit_step(double s);

//Alpha Function
double alpha(double s);

//Resistive Wall Wake Functions
double w_t_rf(double s, double a, double g, double L);
double w_l_rf(double s, double a, double g, double L);

//CSR Wake Function
double w_l_csr(double s, amrex::ParticleReal R, amrex::ParticleReal beam_charge);

//Function to perform FFT-based convolution
void convolve_fft(const std::vector<double>& beam_profile, const std::vector<double>& wake_func, double delta_t, std::vector<double>& result, int padding_factor);

#endif // WAKE_CONVOLUTION_H
