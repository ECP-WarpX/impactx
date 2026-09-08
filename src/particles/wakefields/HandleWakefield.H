/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
#include "ShortRangeWakeElement.H"
#include "WakeConvolution.H"
#include "WakePush.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuParallelReduce.H>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx::particles::wakefields
{
    /** Function to handle wakefield calculations, including charge deposition,
     * mean transverse position calculation, wakefield generation, and convolution.
     *
     * Returns immediately if slice_ds is zero.
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
        if (slice_ds == 0.0) { return; }

        BL_PROFILE("impactx::particles::wakefields::HandleWakefield")

        amrex::ParmParse pp_algo("algo");

        bool csr = false;
        pp_algo.queryAdd("csr", csr);

        bool short_range_wakefields = false;
        pp_algo.queryAdd("short_range_wakefields", short_range_wakefields);

        // Call the CSR bend function
        auto const [element_has_csr, R] = impactx::particles::wakefields::CSRBendElement(element_variant, particle_container.GetRefParticle());

        // Call the function preparing short range wake elements
        auto const short_range_wake = impactx::particles::wakefields::ShortRangeWakeElement(element_variant, particle_container.GetRefParticle());

        // Which effect applies to this very element?
        bool const apply_csr = csr && element_has_csr;
        bool const apply_short_range_wake = short_range_wakefields && short_range_wake.has_wake;

        // Enter loop if wakefields are turned on
        if (apply_csr || apply_short_range_wake)
        {
#ifndef ImpactX_USE_FFT
            throw std::runtime_error(
                std::string(apply_csr ? "algo.csr" : "algo.short_range_wakefields") +
                " was requested but ImpactX was not compiled with FFT support. Recompile with ImpactX_FFT=ON."
            );
#endif

            int csr_bins = 150;
            pp_algo.queryAddWithParser("csr_bins", csr_bins);

            int wake_bins = 150;
            pp_algo.queryAddWithParser("wake_bins", wake_bins);

            // Measure beam size, extract the min, max of particle positions
            [[maybe_unused]] auto const [x_min, y_min, t_min, x_max, y_max, t_max] =
                particle_container.MinAndMaxPositions();

            // Ensure we have a beam with a finite extent in t
            if (t_min == t_max)
            {
                auto const npart = particle_container.TotalNumberOfParticles();
                std::string error_msg = "Coherent Synchrotron Radiation (CSR) "
                                        "modeling was enabled, but ";
                if (npart == 1) {
                    error_msg.append("there is only one particle in the beam. ");
                } else {
                    error_msg.append("the beam is flat in t. ");
                }
                error_msg.append(
                    "CSR modeling is a collective effect and requires "
                    "a particle beam with a finite extent in t. "
                    "Please set the CSR option to false."
                );
                throw std::runtime_error(error_msg);
            }

            // Set parameters for charge deposition
            bool const is_unity_particle_weight = false; // Only true if w = 1
            bool const GetNumberDensity = true;

            int const num_bins = apply_csr ? csr_bins : wake_bins;  // Set resolution
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
                charge_distribution,
                amrex::ParallelDescriptor::IOProcessorNumber(),
                amrex::ParallelDescriptor::Communicator()
            );

            // Pre-size on every rank to the (global) wake length.  On the IO
            // processor this is overwritten by the convolution result of the
            // same length.
            amrex::Gpu::DeviceVector<amrex::Real> convolved_wakefield(num_bins, 0.0);
            if (amrex::ParallelDescriptor::IOProcessor()) {

                // Data for line density (derivative, in the case of CSR) source and 1D wake function
                amrex::Gpu::DeviceVector<amrex::Real> source_data(charge_distribution.size() - 1, 0.0);
                amrex::Gpu::DeviceVector<amrex::Real> wake_function(num_bins*2, 0.0);
                amrex::Real *const dptr_wake_function = wake_function.data();

                // Derivative of the line density is required for the CSR wake
                if (apply_csr){
                    // Call charge density derivative function
                    impactx::particles::wakefields::DerivativeCharge1D(charge_distribution, source_data, bin_size,
                                                                   GetNumberDensity); // Use number derivatives for convolution with CSR
                }

                // Slice mean transverse position is required for the short range structure wake
                if (apply_short_range_wake){
                    // Call the mean transverse position function
                    impactx::particles::wakefields::MeanTransversePosition(particle_container, mean_x, mean_y, bin_min,
                                                                           bin_size, is_unity_particle_weight);
                }

                if (apply_csr){
                    // Construct wake function on 2N support
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
                    convolved_wakefield = impactx::particles::wakefields::convolve_fft(source_data, wake_function, bin_size);
                }
                if (apply_short_range_wake){
                    // TODO:  This section needs to be updated to use the interpolated data from
                    //        short_range_wake.z, short_range_wake.wake_z, etc.
                    // Construct wake function on 2N support
                    auto const wake = short_range_wake;  // trivially copyable view, captured by value
                    amrex::ParallelFor(num_bins*2, [=] AMREX_GPU_DEVICE(int i) {
                        if (i < num_bins) {
                            amrex::Real const s = static_cast<amrex::Real>(i) * bin_size;
                            dptr_wake_function[i] = impactx::particles::wakefields::w_user(s, wake.wake_z, wake.n_data, bin_min, bin_size);
                        }
                        else if (i > num_bins) {
                            amrex::Real const s = static_cast<amrex::Real>(i - 2*num_bins) * bin_size;
                            dptr_wake_function[i] = impactx::particles::wakefields::w_user(s, wake.wake_z, wake.n_data, bin_min, bin_size);
                        }
                    });

                    // Call convolution function
                    convolved_wakefield = impactx::particles::wakefields::convolve_fft(source_data, wake_function, bin_size);
                }

            }

            // Broadcast the global wakefield from the IO processor (which alone
            // computed it) to every MPI rank. Receivers are pre-sized above.
            amrex::ParallelDescriptor::Bcast(
                convolved_wakefield,
                amrex::ParallelDescriptor::IOProcessorNumber(),
                amrex::ParallelDescriptor::Communicator()
            );

            // Check convolution output
            if (print_wakefield && amrex::ParallelDescriptor::IOProcessor())
            {
                // convolved_wakefield is in (non-managed) device memory. Copy it
                // to the host before iterating it on the host for printing.
                amrex::Gpu::HostVector<amrex::Real> h_wakefield(convolved_wakefield.size());
                amrex::Gpu::dtoh_memcpy(h_wakefield.data(), convolved_wakefield.data(),
                                        convolved_wakefield.size() * sizeof(amrex::Real));
                std::cout << "Convolved wakefield: ";
                std::ofstream outfile("convolved_wakefield.txt");
                for (amrex::Real const val : h_wakefield) {
                    std::cout << val << " ";
                    outfile << val << "\n";
                }
                std::cout << "\n";
                outfile.close();
            }

            // Call function to kick particles with wake
            impactx::particles::wakefields::WakePush(particle_container, convolved_wakefield, slice_ds, bin_size, t_min);
        }
    }

} // namespace impactx::particles::wakefields

#endif // HANDLE_WAKEFIELD_H
