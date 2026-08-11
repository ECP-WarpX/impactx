/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
#include <variant>


namespace impactx
{
    void ImpactX::validate ()
    {
        BL_PROFILE("ImpactX::validate");

        // elements
        if (m_lattice->empty())
            throw std::runtime_error("Beamline lattice has zero elements. Not yet initialized?");

        // does a source element at the beginning of the beamline load the
        // beam and, by default, also the reference particle during tracking?
        bool source_loads_beam = false;
        bool source_loads_ref = false;
        if (auto const * source = elements::get_if<elements::Source>(m_lattice->front())) {
            source_loads_beam = true;
            source_loads_ref = source->m_load_ref_particle;
        }

        // reference particle initialized?
        auto const & ref = amr_data->track_particles.m_particle_container->GetRefParticle();
        if (ref.kin_energy_MeV() == 0.0 && !source_loads_ref)
            throw std::runtime_error("The reference particle energy is zero. Not yet initialized?");

        // particles in the beam bunch
        // count particles - if no particles are found in our particle container, then a lot of
        // AMReX routines over ParIter won't work, and we have nothing to do here anyway
        {
            int const nLevelPC = amr_data->finestLevel();
            amrex::Long nParticles = 0;
            for (int lev = 0; lev <= nLevelPC; ++lev) {
                nParticles += amr_data->track_particles.m_particle_container->NumberOfParticlesAtLevel(lev);
            }
            if (nParticles == 0 && !source_loads_beam)
            {
                throw std::runtime_error(
                    "No particles found. "
                    "Cannot track particles without an initialized beam. "
                    "Did you forget to call sim.add_particles ?"
                );
            }
        }
    }
} // namespace impactx
