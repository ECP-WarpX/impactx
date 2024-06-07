#ifndef WAKE_CONVOLUTION_H
#define WAKE_CONVOLUTION_H

#include "../ImpactXParticleContainer.H"
#include "../../initialization/InitDistribution.H"

using namespace amrex;

//Constants
constexpr double Z0 = 377;  //Free space impedance [Ohm]
constexpr double alpha_1 = 0.4648;

//Step Function
double unit_step(double s);

//Alpha Function
double alpha(double s);

//Resistive Wall Wake Functions
double W_T_RF(double s, double a, double g, double L);
double W_L_RF(double s, double a, double g, double L);

//CSR Wake Function
double W_L_CSR(double s, amrex::ParticleReal R, amrex::ParticleReal beam_charge);

//Function to perform FFT-based convolution
void Convolve_FFT(const std::vector<double>& beam_profile, const std::vector<double>& wake_func, double delta_t, std::vector<double>& result, int padding_factor);

#endif // WAKE_CONVOLUTION_H
