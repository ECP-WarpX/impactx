/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef CHARGE_BINNING_H
#define CHARGE_BINNING_H

#include "../ImpactXParticleContainer.H"

using namespace amrex;

//Function to calculate charge density profile
void DepositCharge1D(impactx::ImpactXParticleContainer& myspc, Real* dptr_data, int num_bins, Real bin_min, Real bin_size, bool is_unity_particle_weight);

//Function to calculate number, charge density slopes
void DerivativeCharge1D(double* charge_distribution, double* slopes, int num_bins, Real bin_size, bool GetNumberDensity);

#endif // CHARGE_BINNING_H
