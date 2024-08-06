/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
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
double w_l_csr(double s, amrex::ParticleReal R);

//Function to perform FFT-based convolution
void convolve_fft(double* beam_profile, double* wake_func, int beam_profile_size, int wake_func_size, double delta_t, double* result, int padding_factor);

#endif // WAKE_CONVOLUTION_H
