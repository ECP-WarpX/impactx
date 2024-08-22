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
#include "particles/ImpactXParticleContainer.H"

#ifdef ImpactX_USE_FFT
#include <ablastr/math/fft/AnyFFT.H>
#endif

#include <AMReX_REAL.H>

#include <algorithm>
#include <cmath>
#ifndef ImpactX_USE_FFT
#include <stdexcept>
#endif

namespace impactx::particles::wakefields
{
    amrex::Real alpha (amrex::Real s)
    {
        using namespace amrex::literals;

        return 1_rt - impactx::particles::wakefields::alpha_1 * std::sqrt(s) - (1_rt - 2_rt * impactx::particles::wakefields::alpha_1) * s;
    }

    amrex::Real w_t_rf (
        amrex::Real s,
        amrex::Real a,
        amrex::Real g,
        amrex::Real L
    )
    {
        using namespace amrex::literals;

        amrex::Real const s0 = (0.169_rt * std::pow(a, 1.79_rt) * std::pow(g, 0.38_rt)) / std::pow(L, 1.17_rt);
        amrex::Real const term = std::sqrt(std::abs(s) / s0) * std::exp(-std::sqrt(std::abs(s) / s0));
        return (4 * impactx::particles::wakefields::Z0 * ablastr::constant::SI::c * s0 * unit_step(s)) / (amrex::Real(M_PI) * std::pow(a, 4)) * term;
    }

    amrex::Real w_l_rf (
        amrex::Real s,
        amrex::Real a,
        amrex::Real g,
        amrex::Real L
    )
    {
        using namespace amrex::literals;

        amrex::Real const s00 = g * std::pow((a / (alpha(g / L) * L)), 2) / 8.0_rt;
        return (impactx::particles::wakefields::Z0 * ablastr::constant::SI::c * unit_step(s) * std::exp(-std::sqrt(std::abs(s) / s00))) / (amrex::Real(M_PI) * std::pow(a, 2));
    }

    amrex::Gpu::DeviceVector<amrex::Real>
    convolve_fft (
        amrex::Gpu::DeviceVector<amrex::Real> const & beam_profile_slope,
        amrex::Gpu::DeviceVector<amrex::Real> const & wake_func,
        amrex::Real delta_t,
        int padding_factor
    )
    {
    #ifdef ImpactX_USE_FFT
        int const beam_profile_slope_size = beam_profile_slope.size();
        int const wake_func_size = wake_func.size();

        // Length of convolution result (complex) is N/2+1
        // e.g., https://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html
        int const signal_size = padding_factor * std::max(beam_profile_slope_size, wake_func_size);  // zero-pad slope by one element
        int const complex_size = signal_size / 2 + 1;

        // Allocate memory for FFT inputs and outputs
        using ablastr::math::anyfft::Complex;

        // Zero-pad the input (signal) arrays to equal length
        amrex::Gpu::DeviceVector<amrex::Real> in1(signal_size, 0.0);
        amrex::Gpu::DeviceVector<amrex::Real> in2(signal_size, 0.0);
        amrex::Real * const dptr_in1 = in1.data();
        amrex::Real * const dptr_in2 = in2.data();
        amrex::Real const * const dptr_beam_profile_slope = beam_profile_slope.data();
        amrex::Real const * const dptr_wake_func = wake_func.data();
        amrex::ParallelFor(signal_size, [=] AMREX_GPU_DEVICE(int i)
        {
            if (i < beam_profile_slope_size)
            {
                dptr_in1[i] = dptr_beam_profile_slope[i];
            }
            else
            {
                dptr_in1[i] = 0;
            }

            if (i < wake_func_size)
            {
                dptr_in2[i] = dptr_wake_func[i];
            }
            else
            {
                dptr_in2[i] = 0;
            }
        });

        // Define Forward FFT
        amrex::Gpu::DeviceVector<Complex> out1(complex_size);
        amrex::Gpu::DeviceVector<Complex> out2(complex_size);
        // TODO: n does not change usually, so we can keep the plans alive over the simulation
        //       runtime. To do that, we can make this function a functor class.
        auto p1 = ablastr::math::anyfft::CreatePlan(
            amrex::IntVect{signal_size}, in1.data(), out1.data(), ablastr::math::anyfft::direction::R2C, 1

        );
        auto p2 = ablastr::math::anyfft::CreatePlan(
            amrex::IntVect{signal_size}, in2.data(), out2.data(), ablastr::math::anyfft::direction::R2C, 1
        );

        // Perform Forward FFT - Convert inputs into frequency domain
        // Gives out1,out2, the FFT-transformed input arrays of in1, in2
        ablastr::math::anyfft::Execute(p1);
        ablastr::math::anyfft::Execute(p2);

        // Perform FFT Multiplication - FFT Element-wise multiplication in frequency space
        amrex::Gpu::DeviceVector<Complex> conv_result(complex_size);
        Complex * const dptr_conv_result = conv_result.data();
        Complex const * const dptr_out1 = out1.data();
        Complex const * const dptr_out2 = out2.data();
        amrex::ParallelFor(complex_size, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            using ablastr::math::anyfft::multiply;
            multiply(dptr_conv_result[i], dptr_out1[i], dptr_out2[i]);
        });

        // Define Backward FFT - Revert from frequency domain to time/space domain
        amrex::Gpu::DeviceVector<amrex::Real> result(signal_size, 0.0);
        // TODO: n does not change usually, so we can keep the plans alive over the simulation
        //       runtime. To do that, we can make this function a functor class.
        amrex::Real * const dptr_result = result.data();
        auto p3 = ablastr::math::anyfft::CreatePlan(
            amrex::IntVect{signal_size}, dptr_result, dptr_conv_result, ablastr::math::anyfft::direction::C2R, 1
        );

        // Perform Backward FFT
        ablastr::math::anyfft::Execute(p3);

        // Normalize result by the output size and multiply result by bin size
        amrex::ParallelFor(signal_size, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            dptr_result[i] = dptr_result[i] / signal_size * delta_t;
        });

        // Clean up intermediate declarations
        // TODO: n does not change usually, so we can keep the plans alive over the simulation
        //       runtime. To do that, we can make this function a functor class.
        ablastr::math::anyfft::DestroyPlan(p1);
        ablastr::math::anyfft::DestroyPlan(p2);
        ablastr::math::anyfft::DestroyPlan(p3);

        return result;
    #else
        throw std::runtime_error("convolve_fft: To use this function, recompile with ImpactX_FFT=ON.");

        return amrex::Gpu::DeviceVector<amrex::Real>();
    #endif
    }
}
