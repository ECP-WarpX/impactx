#ifndef CHARGE_BINNING_H
#define CHARGE_BINNING_H

#include "../ImpactXParticleContainer.H"

using namespace amrex;

//Function to calculate charge density profile
void DepositCharge1D(impactx::ImpactXParticleContainer& myspc, Real* dptr_data, int num_bins, Real bin_min, Real bin_size, bool is_unity_particle_weight);

//Function to calculate number, charge density slopes
void DerivativeCharge1D(std::vector<double>& charge_distribution, std::vector<double>& slopes, int num_bins, Real bin_size, bool GetNumberDensity);

#endif // CHARGE_BINNING_H
