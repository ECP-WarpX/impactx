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

#include <AMReX_REAL.H>


namespace impactx::particles::wakefields
{
    // Constants
    constexpr amrex::Real Z0 = 377;  //!< Free space impedance [Ohm]

    /** Wake function constant [unitless]
     *
     * Used in eq. (21) of:
     * K. L.F. Bane, "Short-Range Dipole Wakefields in Accelerating Structures for the NLC," SLAC-PUB-9663, 2003
     */
    constexpr amrex::Real alpha_1 = 0.4648;

    /** Step Function (Heavy-Side)
     *
     * @param s value along s [m]
     * @return returns 1 if s>0, else 0
     */
    amrex::Real unit_step (amrex::Real s);

    /** Alpha Function
     *
     * Eq. (21) in:
     * K. L.F. Bane, "Short-Range Dipole Wakefields in Accelerating Structures for the NLC," SLAC-PUB-9663, 2003
     *
     * @param s value along s [m]
     * @return value of alpha(s) in Eq. (21) of Bane, 2003
     */
    amrex::Real alpha (amrex::Real s);

    /** Resistive Wall Wake Function (transverse)
     *
     * Eq. (17) in:
     * K. L.F. Bane, "Short-Range Dipole Wakefields in Accelerating Structures for the NLC," SLAC-PUB-9663, 2003
     *
     * @param[in] s value along s [m]
     * @param[in] a iris radius [m]
     * @param[in] g gap [m]
     * @param[in] L period length [m]
     * @return transverse wake field [V/m/pC/mm]
     */
    amrex::Real w_t_rf (
        amrex::Real s,
        amrex::Real a,
        amrex::Real g,
        amrex::Real L
    );

    /** Resistive Wall Wake Function (longitudinal)
     *
     * Eq. (18) in:
     * K. L.F. Bane, "Short-Range Dipole Wakefields in Accelerating Structures for the NLC," SLAC-PUB-9663, 2003
     *
     * @param[in] s value along s [m]
     * @param[in] a iris radius [m]
     * @param[in] g gap [m]
     * @param[in] L period length [m]
     * @return longitudinal wake field [V/m/pC/mm]
     */
    amrex::Real w_l_rf (
        amrex::Real s,
        amrex::Real a,
        amrex::Real g,
        amrex::Real L
    );

    /** CSR Wake Function
     *
     * Eq. (28) in:
     * E.L. Saldin et al., "On the coherent radiation of an electron bunch moving in an arc of a circle",
     * NiMA Volume 398, Issues 2–3, Pages 373-394, 21 October (1997)
     * https://doi.org/10.1016/S0168-9002(97)00822-X
     *
     * @param[in] s value along s [m]
     * @param[in] R band radius [m]
     * @return wake functions in [V*pc/mm]
     */
    amrex::Real w_l_csr (
        amrex::Real s,
        amrex::Real R
    );

    /** Perform FFT-based convolution
     *
     * @param[in] beam_profile_slope number density slope along s [1/m]
     * @param[in] wake_func wake functions in along s [V*pc/mm]
     * @param[in] beam_profile_size number of bins in beam_profile_slope
     * @param[in] wake_func_size number of bins in wake_func
     * @param[in] delta_t size of a bin in wake_func [m]
     * @param[out] result FFT convolution of beam_profile_slope & wake_func
     * @param[in] padding_factor N times zero-padding of the input arrays' lengths
     */
    void convolve_fft (
        amrex::Real const * beam_profile_slope,
        amrex::Real const * wake_func,
        int beam_profile_size,
        int wake_func_size,
        amrex::Real delta_t,
        amrex::Real* result,
        int padding_factor
    );

} // namespace impactx::particles::wakefields

#endif // WAKE_CONVOLUTION_H
