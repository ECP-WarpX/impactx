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
#ifndef ImpactX_USE_FFT
#include <stdexcept>
#endif
#include <vector>


namespace impactx::particles::wakefields
{
    /** Function to handle CSR bending process including charge deposition,
     * mean transverse position calculation, wakefield generation, and convolution.
     *
     * @param[in] particle_container the particle species container
     * @param[in] element_variant variant type of the lattice element
     * @param[in] slice_ds slice spacing along s
     * @param[in] print_wakefield for debugging: print the wakefield to convolved_wakefield.txt
     */
    template <typename T_Element>
    void HandleWakefield (
        impactx::ImpactXParticleContainer& particle_container,
        T_Element const& element_variant,
        amrex::Real slice_ds,
        bool print_wakefield = false
    )
    {
        BL_PROFILE("impactx::particles::wakefields::HandleWakefield")

        amrex::ParmParse pp_algo("algo");
        bool csr = false;
        pp_algo.queryAdd("csr", csr);

        // Call the CSR bend function
        auto const [element_has_csr, R] = impactx::particles::wakefields::CSRBendElement(element_variant, particle_container.GetRefParticle());

        // Enter loop if lattice has bend element
        if (csr && element_has_csr)
        {
#ifndef ImpactX_USE_FFT
            throw std::runtime_error("algo.csr was requested but ImpactX was not compiled with FFT support. Recompile with ImpactX_FFT=ON.");
#endif

            int csr_bins = 150;
            pp_algo.queryAdd("csr_bins", csr_bins);

            // Measure beam size, extract the min, max of particle positions
            [[maybe_unused]] auto const [x_min, y_min, t_min, x_max, y_max, t_max] =
                particle_container.MinAndMaxPositions();

            // Set parameters for charge deposition
            bool const is_unity_particle_weight = false; // Only true if w = 1
            bool const GetNumberDensity = true;

            int const num_bins = csr_bins;  // Set resolution
            amrex::Real const bin_min = t_min;
            amrex::Real const bin_max = t_max;
            amrex::Real const bin_size = (bin_max - bin_min) / (num_bins - 1);  // number of evaluation points

            // Allocate memory for the charge profile
            amrex::Gpu::DeviceVector<amrex::Real> charge_distribution(num_bins + 1, 0.0);
            amrex::Gpu::DeviceVector<amrex::Real> mean_x(num_bins, 0.0);
            amrex::Gpu::DeviceVector<amrex::Real> mean_y(num_bins, 0.0);

            // Call charge deposition function
            impactx::particles::wakefields::DepositCharge1D(particle_container, charge_distribution, bin_min, bin_size, is_unity_particle_weight);

            // Sum up all partial charge histograms to one MPI process to calculate the global wakefield.
            // Once calculated, we will distribute convolved_wakefield back to every MPI process.
            amrex::ParallelReduce::Sum(
                charge_distribution.data(),
                charge_distribution.size(),
                amrex::ParallelDescriptor::IOProcessorNumber(),
                amrex::ParallelDescriptor::Communicator()
            );

            amrex::Gpu::DeviceVector<amrex::Real> convolved_wakefield;
            if (amrex::ParallelDescriptor::IOProcessor()) {
                // Call the mean transverse position function
                impactx::particles::wakefields::MeanTransversePosition(particle_container, mean_x, mean_y, bin_min,
                                                                       bin_size, is_unity_particle_weight);

                // Call charge density derivative function
                amrex::Gpu::DeviceVector<amrex::Real> slopes(charge_distribution.size() - 1, 0.0);
                impactx::particles::wakefields::DerivativeCharge1D(charge_distribution, slopes, bin_size,
                                                                   GetNumberDensity); // Use number derivatives for convolution with CSR

                // Construct CSR wake function on 2N support
                amrex::Gpu::DeviceVector<amrex::Real> wake_function(num_bins*2, 0.0);
                amrex::Real *const dptr_wake_function = wake_function.data();
                auto const dR = R;  // for NVCC capture
                amrex::ParallelFor(num_bins*2, [=] AMREX_GPU_DEVICE(int i) {
                    if (i < num_bins) {
                        amrex::Real const s = static_cast<amrex::Real>(i) * bin_size;
                        dptr_wake_function[i] = impactx::particles::wakefields::w_l_csr(s, dR, bin_size);
                    }
                    else if (i > num_bins) {
                        amrex::Real const s = static_cast<amrex::Real>(i - 2*num_bins) * bin_size;
                        dptr_wake_function[i] = impactx::particles::wakefields::w_l_csr(s, dR, bin_size);
                    }
                });

                // Call convolution function
                convolved_wakefield = impactx::particles::wakefields::convolve_fft(slopes, wake_function, bin_size);
            }

            // Broadcast the global wakefield to every MPI rank
            amrex::ParallelDescriptor::Bcast(
                convolved_wakefield.data(),
                convolved_wakefield.size(),
                amrex::ParallelDescriptor::IOProcessorNumber()
            );

            // Check convolution output
            if (print_wakefield && amrex::ParallelDescriptor::IOProcessor())
            {
                std::cout << "Convolved wakefield: ";
                std::ofstream outfile("convolved_wakefield.txt");
                for (double const i : convolved_wakefield) {
                    std::cout << i << " ";
                    outfile << i << std::endl;
                }
                std::cout << std::endl;
                outfile.close();
            }

            // Call function to kick particles with wake
            impactx::particles::wakefields::WakePush(particle_container, convolved_wakefield, slice_ds, bin_size, t_min);
        }
    }

} // namespace impactx::particles::wakefields

#endif // HANDLE_WAKEFIELD_H
