/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef HANDLE_WAKEFIELD_H
#define HANDLE_WAKEFIELD_H

#include "particles/ImpactXParticleContainer.H"
#include "ChargeBinning.H"
#include "CSRBendElement.H"
#include "WakeConvolution.H"
#include "WakePush.H"

#include <cmath>
#include <vector>

using namespace amrex;

namespace impactx::particles::wakefields
{

    /** Function to handle CSR bending process including charge deposition,
     * mean transverse position calculation, wakefield generation, and convolution.
     *
     * @param[in] particle_container the particle species container
     * @param[in] R radius of curvature of the bend element
     * @param[in] element_has_csr boolean indicating if the element has CSR
     * @param[in] element_variant variant type of the lattice element
     * @param[in] csr boolean indicating if CSR effects are included
     * @param[in] csr_bins number of bins (resolution) for CSR calculations
     * @param[in] slice_ds slice spacing along s
     */

    template <typename T_Element>
    void HandleWakefield(
        impactx::ImpactXParticleContainer* particle_container,
        amrex::Real R,
        bool element_has_csr,
        T_Element const& element_variant,
        bool csr,
        int csr_bins,
        amrex::Real slice_ds)
    {

        // Call the CSR bend function
        impactx::particles::wakefields::CSRBendElement(R, element_has_csr, element_variant);
        //impactx::particles::wakefields::HandleElementVariant(R, element_has_csr, element_variant);

        std::cout << "My boolean for element has csr is:" << element_has_csr << std::endl;

        // Enter loop if lattice has bend element
        if (csr && element_has_csr)
        {
            // Measure beam size, extract the min, max of particle positions
            auto const [x_min, y_min, t_min, x_max, y_max, t_max] =
                particle_container->MinAndMaxPositions();

            using amrex::Real;

            // Set parameters for charge deposition
            bool is_unity_particle_weight = false; // Only true if w = 1
            bool GetNumberDensity = true;

            int padding_factor = 1; // Set amount of zero-padding
            int num_bins = csr_bins;  // Set resolution
            Real bin_min = t_min;
            Real bin_max = t_max;
            Real bin_size = (bin_max - bin_min) / num_bins;

            // Allocate memory for the charge profile
            Real* dptr_data = new Real[num_bins]();
            Real* mean_x = new Real[num_bins]();
            Real* mean_y = new Real[num_bins]();

            // Call charge deposition function
            impactx::particles::wakefields::DepositCharge1D(*particle_container, dptr_data, num_bins, bin_min, bin_size, is_unity_particle_weight);

            // Call the mean transverse position function
            impactx::particles::wakefields::MeanTransversePosition(*particle_container, mean_x, mean_y, num_bins, bin_min, bin_size, is_unity_particle_weight);

            // Call charge density derivative function
            std::vector<amrex::Real> charge_distribution(dptr_data, dptr_data + num_bins);
            std::vector<amrex::Real> slopes(num_bins - 1);
            impactx::particles::wakefields::DerivativeCharge1D(charge_distribution.data(), slopes.data(), num_bins, bin_size, GetNumberDensity); // Use number derivatives for convolution with CSR

            // Call CSR wake function
            std::vector<amrex::Real> wake_function(num_bins);
            for (int i = 0; i < num_bins; ++i)
            {
                amrex::Real s = bin_min + i * bin_size;
                wake_function[i] = impactx::particles::wakefields::w_l_csr(s, R);
            }

            // Call convolution function
            std::vector<amrex::Real> convoluted_wakefield(padding_factor * (2 * num_bins - 1));
            impactx::particles::wakefields::convolve_fft(slopes.data(), wake_function.data(), slopes.size(), wake_function.size(), bin_size, convoluted_wakefield.data(), padding_factor);

            // Check convolution output
            std::cout << "Convoluted wakefield: ";
            std::ofstream outfile("convoluted_wakefield.txt");
            for (int i = 0; i < int(convoluted_wakefield.size()); ++i)
            {
                std::cout << convoluted_wakefield[i] << " ";
                outfile << convoluted_wakefield[i] << std::endl;
            }
            std::cout << std::endl;
            outfile.close();
            delete[] dptr_data;

            // Call function to kick particles with wake
            impactx::particles::wakefields::WakePush(*particle_container, convoluted_wakefield, slice_ds, bin_size, t_min, padding_factor);
        }
    }

} // namespace impactx::particles::wakefields

#endif // HANDLE_WAKEFIELD_H
