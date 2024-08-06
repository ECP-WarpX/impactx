/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "WakeConvolution.H"
#include "particles/ImpactXParticleContainer.H" //Includes all necessary AMReX headers
#include "initialization/InitDistribution.H"

#ifdef ImpactX_USE_FFT
#include <ablastr/math/fft/AnyFFT.H>
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

/*
Wake Functions:

Define the longitudinal and transverse resistive wall wake functions,
the longitudinal CSR wake function, and their associated functions

*/

//Step Function
double unit_step(double s)
{
    return s >= 0 ? 1.0 : 0.0; //If true return 1, else 0
}

//Alpha Function
double alpha(double s)
{
    return 1.0 - alpha_1 * std::sqrt(s) - (1.0 - 2 * alpha_1) * s;
}

//Resistive Wall Wake Functions

double w_t_rf(double s, double a, double g, double L)
{
    double s0 = (0.169 * std::pow(a, 1.79) * std::pow(g, 0.38)) / std::pow(L, 1.17);
    double term = std::sqrt(std::abs(s) / s0) * std::exp(-std::sqrt(std::abs(s) / s0));
    return (4 * Z0 * ablastr::constant::SI::c * s0 * unit_step(s)) / (M_PI * std::pow(a, 4)) * term;
}

double w_l_rf(double s, double a, double g, double L)
{
    double s00 = g * std::pow((a / (alpha(g / L) * L)), 2) / 8.0;
    return (Z0 * ablastr::constant::SI::c * unit_step(s) * std::exp(-std::sqrt(std::abs(s) / s00))) / (M_PI * std::pow(a, 2));
}

//CSR Wake Function

double w_l_csr(double s, amrex::ParticleReal R)
{
    double rc = std::pow(ablastr::constant::SI::q_e, 2) / (4 * M_PI * ablastr::constant::SI::ep0 * ablastr::constant::SI::m_e * std::pow(ablastr::constant::SI::c, 2));
    double kappa = (2 * rc * ablastr::constant::SI::m_e * std::pow(ablastr::constant::SI::c, 2)) / std::pow(3, 1.0/3.0) / std::pow(R, 2.0/3.0);

    return - (kappa * unit_step(s)) / std::pow(std::abs(s), 1.0/3.0);
}

//Convolution Function

void convolve_fft(double* beam_profile, double* wake_func, int beam_profile_size, int wake_func_size, double delta_t, double* result, int padding_factor)
{
#ifdef ImpactX_USE_FFT
    //Length of convolution result
    int original_n = beam_profile_size + wake_func_size - 1; //Output size is n = 2N - 1, where N = size of signals 1,2

    //Add padding factor to control amount of zero-padding
    int n = static_cast<int>(original_n * padding_factor);

    //Allocate memory for FFT inputs and outputs
    using ablastr::math::anyfft::Complex;
    double *in1 = (double*) malloc(sizeof(double) * n); //Allocate memory for 'n' real numbers for inputs and complex outputs
    double *in2 = (double*) malloc(sizeof(double) * n);
    Complex *out1 = (Complex*) malloc(sizeof(Complex) * n);
    Complex *out2 = (Complex*) malloc(sizeof(Complex) * n);
    Complex *conv_result = (Complex*) malloc(sizeof(Complex) * n);
    double *out3 = (double*) malloc(sizeof(double) * n);

    //Zero-pad the input arrays to be the size of the convolution output length 'n'
    for (int i = 0; i < n; ++i)
    {
        if (i < beam_profile_size)
        {
            in1[i] = std::isfinite(beam_profile[i]) ? beam_profile[i] : 0.0; //Print NaN was produced if 0
        }
        else
        {
            in1[i] = 0.0;
        }

        if (i < wake_func_size)
        {
            in2[i] = std::isfinite(wake_func[i]) ? wake_func[i] : 0.0; //Print NaN was produced if 0
        }
        else
        {
            in2[i] = 0.0;
        }
    }

    // Define Forward FFT
    auto p1 = ablastr::math::anyfft::CreatePlan(
        amrex::IntVect{n}, in1, out1, ablastr::math::anyfft::direction::R2C, 1
    );
    auto p2 = ablastr::math::anyfft::CreatePlan(
        amrex::IntVect{n}, in2, out2, ablastr::math::anyfft::direction::R2C, 1
    );

    // Perform Forward FFT - Convert inputs into frequency domain
    // Gives out1,out2, the FFT-transformed input arrays of in1, in2
    ablastr::math::anyfft::Execute(p1);
    ablastr::math::anyfft::Execute(p2);

    //Perform FFT Multiplication - FFT Element-wise multiplication in frequency space
    for (int i = 0; i < n; ++i)
    {
        /*
        The multiplication of 2 complex numbers is given by:

        A * B = (a_real + i * a_imag)(b_real + i * b_imag)
        A * B = (a_real * b_real - a_imag * b_imag) + i * (a_real * b_imag + a_imag * b_real)

        where out[i][0] is the real part and out[i][1] is the imaginary part of out1, out2
        */

        conv_result[i][0] = out1[i][0] * out2[i][0] - out1[i][1] * out2[i][1]; //Real part of convolution
        conv_result[i][1] = out1[i][0] * out2[i][1] + out1[i][1] * out2[i][0]; //Imaginary part of convolution
    }

    // Define Backward FFT - Revert from frequency domain to time/space domain
    auto p3 = ablastr::math::anyfft::CreatePlan(
        amrex::IntVect{n}, out3, conv_result, ablastr::math::anyfft::direction::C2R, 1
    );

    // Perform Backward FFT
    ablastr::math::anyfft::Execute(p3);

    //Normalize result by the output size and multiply result by bin size
    for (int i = 0; i < n; ++i)
    {
        result[i] = out3[i] / n * delta_t;
    }

    //Clean up intermediate declarations
    ablastr::math::anyfft::DestroyPlan(p1);
    ablastr::math::anyfft::DestroyPlan(p2);
    ablastr::math::anyfft::DestroyPlan(p3);
    free(in1);
    free(in2);
    free(out1);
    free(out2);
    free(out3);
    free(conv_result);
#else
    throw std::runtime_error("convolve_fft: To use this function, recompile with ImpactX_FFT=ON.");
#endif
}
