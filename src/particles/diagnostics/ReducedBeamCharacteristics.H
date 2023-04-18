/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_REDUCED_BEAM_CHARACTERISTICS_H
#define IMPACTX_REDUCED_BEAM_CHARACTERISTICS_H

#include "particles/ImpactXParticleContainer.H"
#include "ImpactX.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

using namespace amrex;

namespace impactx::diagnostics
{
    /** Struct containing reduced beam characteristics data
     */
     struct ReducedBeamCharacteristics{
         ParticleReal s;  // TODO: remove when the output gets rerouted to openPMD
         ParticleReal beta_gamma;  // TODO: remove when the output gets rerouted to openPMD
         ParticleReal x_mean;
         ParticleReal y_mean;
         ParticleReal t_mean;
         ParticleReal sig_x;
         ParticleReal sig_y;
         ParticleReal sig_t;
         ParticleReal ux_mean;
         ParticleReal uy_mean;
         ParticleReal ut_mean;
         ParticleReal sig_px;
         ParticleReal sig_py;
         ParticleReal sig_pt;
         ParticleReal emittance_x;
         ParticleReal emittance_y;
         ParticleReal emittance_t;
         ParticleReal alpha_x;
         ParticleReal alpha_y;
         ParticleReal alpha_t;
         ParticleReal beta_x;
         ParticleReal beta_y;
         ParticleReal beta_t;
         ParticleReal charge;  // TODO: remove when the output gets rerouted to openPMD
     };

    /** Compute momenta of the beam distribution and print them out
      */
    ReducedBeamCharacteristics compute_reduced_beam_characteristics (ImpactXParticleContainer const & pc);

} // namespace impactx

#endif // IMPACTX_REDUCED_BEAM_CHARACTERISTICS
