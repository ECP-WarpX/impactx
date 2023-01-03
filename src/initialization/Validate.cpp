/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_INT.H>

#include <stdexcept>


namespace impactx
{
    void ImpactX::validate ()
    {
        BL_PROFILE("ImpactX::validate");

        // reference particle initialized?
        auto const & ref = m_particle_container->GetRefParticle();
        if (ref.energy_MeV() == 0.0)
            throw std::runtime_error("The reference particle energy is zero. Not yet initialized?");

        // particles in the beam bunch
        // count particles - if no particles are found in our particle container, then a lot of
        // AMReX routines over ParIter won't work, and we have nothing to do here anyway
        {
            int const nLevelPC = finestLevel();
            amrex::Long nParticles = 0;
            for (int lev = 0; lev <= nLevelPC; ++lev) {
                nParticles += m_particle_container->NumberOfParticlesAtLevel(lev);
            }
            if (nParticles == 0)
                throw std::runtime_error("No particles found. Cannot run evolve without a beam.");
            if (nParticles == 1)
                throw std::runtime_error("Only one particle found. This is not yet supported: https://github.com/ECP-WarpX/impactx/issues/44");
        }

        // elements
        if (m_lattice.size() == 0u)
            throw std::runtime_error("Beamline lattice has zero elements. Not yet initialized?");
    }
} // namespace impactx
