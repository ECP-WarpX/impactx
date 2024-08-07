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

namespace impactx::particles::wakefields
{
    //Constants
    constexpr amrex::Real Z0 = 377;  //Free space impedance [Ohm]
    constexpr amrex::Real alpha_1 = 0.4648;

    //Step Function
    amrex::Real unit_step(amrex::Real s);

    //Alpha Function
    amrex::Real alpha(amrex::Real s);

    //Resistive Wall Wake Functions
    amrex::Real w_t_rf(amrex::Real s, amrex::Real a, amrex::Real g, amrex::Real L);
    amrex::Real w_l_rf(amrex::Real s, amrex::Real a, amrex::Real g, amrex::Real L);

    //CSR Wake Function
    amrex::Real w_l_csr(amrex::Real s, amrex::ParticleReal R);

    //Function to perform FFT-based convolution
    void convolve_fft(amrex::Real* beam_profile, amrex::Real* wake_func, int beam_profile_size, int wake_func_size, amrex::Real delta_t, amrex::Real* result, int padding_factor);
} //namespace impactx::particles::wakefields

#endif //WAKE_CONVOLUTION_H
