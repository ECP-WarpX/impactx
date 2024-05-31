#ifndef WAKE_CONVOLUTION_H
#define WAKE_CONVOLUTION_H

#include "ImpactXParticleContainer.H"

using namespace amrex;

//Constants
constexpr double e = 1.602e-19;  // Elementary charge [C]
constexpr double c = 3.0e8;  //Speed of light [m/s]
constexpr double m_e = 9.109e-31;  //Electron mass [kg]
constexpr double Z0 = 377;  //Free space impedance [Ohm]
constexpr double epsilon_0 = 8.854e-12;  //Vacuum permittivity [F/m]
constexpr double alpha_1 = 0.4648;

//Step Function
double unit_step(double s);

//Alpha Function
double alpha(double s);

//Resistive Wall Wake Functions
double W_T_RF(double s, double a, double g, double L);
double W_L_RF(double s, double a, double g, double L);

//CSR Wake Function
double W_L_CSR(double s, double R);

//Function to perform FFT-based convolution
void Convolve_FFT(const std::vector<double>& beam_profile, const std::vector<double>& wake_func, double delta_t, std::vector<double>& result);

#endif // WAKE_CONVOLUTION_H