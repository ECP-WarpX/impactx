#include "WakeConvolution.H"
#include "../ImpactXParticleContainer.H" //Includes all necessary AMReX headers
#include "../../initialization/InitDistribution.H"
 
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fftw3.h> //Fastest Fourier Transform in the West

using namespace amrex;

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

double W_T_RF(double s, double a, double g, double L)
{
    double s0 = (0.169 * std::pow(a, 1.79) * std::pow(g, 0.38)) / std::pow(L, 1.17);
    double term = std::sqrt(std::abs(s) / s0) * std::exp(-std::sqrt(std::abs(s) / s0));
    return (4 * Z0 * c * s0 * unit_step(s)) / (M_PI * std::pow(a, 4)) * term;
}

double W_L_RF(double s, double a, double g, double L)
{
    double s00 = g * std::pow((a / (alpha(g / L) * L)), 2) / 8.0;
    return (Z0 * c * unit_step(s) * std::exp(-std::sqrt(std::abs(s) / s00))) / (M_PI * std::pow(a, 2));
}

//CSR Wake Function

double W_L_CSR(double s, amrex::ParticleReal R, amrex::ParticleReal beam_charge)
{
    double N = beam_charge / e;
    double rc = e * e / (4 * M_PI * epsilon_0 * m_e * std::pow(c, 2));
    double kappa = (2 * rc * m_e * std::pow(c, 2)) / std::pow(3, 1.0/3.0) / std::pow(R, 2.0/3.0);

    return - (N * kappa * unit_step(s)) / std::pow(std::abs(s), 1.0/3.0);
}

//Convolution Function

void Convolve_FFT(const std::vector<double>& beam_profile, const std::vector<double>& wake_func, double delta_t, std::vector<double>& result, int padding_factor)
{
    //Length of convolution result
    int original_n = beam_profile.size() + wake_func.size() - 1; //Output size is n = 2N - 1, where N = size of signals 1,2

    //Add padding factor to control amount of zero-padding
    int n = static_cast<int>(original_n * padding_factor);

    //Allocate memory for FFTW inputs and outputs
    fftw_complex *in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n); //Allocate memory for 'n' complex numbers and
    fftw_complex *in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n); //Convert void pointer --> fftw_complex pointer
    fftw_complex *out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex *out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex *conv_result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    //Zero-pad the input arrays to be the size of the convolution output length 'n'
    for (int i = 0; i < n; ++i)
    {
        if (i < beam_profile.size())
        {
            in1[i][0] = std::isfinite(beam_profile[i]) ? beam_profile[i] : 0.0; //Ensure no nans
            in1[i][1] = 0.0;
        }
        else
        {
            in1[i][0] = 0.0;
            in1[i][1] = 0.0;
        }

        if (i < wake_func.size())
        {
            in2[i][0] = std::isfinite(wake_func[i]) ? wake_func[i] : 0.0; //Ensure no nans
            in2[i][1] = 0.0;
        }
        else
        {
            in2[i][0] = 0.0;
            in2[i][1] = 0.0;
        }
    }

    //Define Forward FFT
    fftw_plan p1 = fftw_plan_dft_1d(n, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE); //sign = forward or backward; flag = estimate or measure
    fftw_plan p2 = fftw_plan_dft_1d(n, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

    //Perform Forward FFT - Convert inputs into frequency domain
    fftw_execute(p1); //Gives out1,out2, the FFT-transformed input arrays of in1, in2
    fftw_execute(p2);

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

    //Define Backward FFT - Revert from frequency domain to time/space domain
    fftw_plan p3 = fftw_plan_dft_1d(n, conv_result, out1, FFTW_BACKWARD, FFTW_ESTIMATE);

    //Perform Backward FFT
    fftw_execute(p3);

    //Normalize result by the output size and multiply result by bin size
    result.resize(n);
    for (int i = 0; i < n; ++i)
    {
        result[i] = out1[i][0] / n * delta_t;
    }

    //Clean up intermediate declarations
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_free(in1);
    fftw_free(in2);
    fftw_free(out1);
    fftw_free(out2);
    fftw_free(conv_result);
}
