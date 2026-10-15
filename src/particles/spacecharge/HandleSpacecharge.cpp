/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "HandleSpacecharge.H"

#include "initialization/Algorithms.H"
#include "initialization/AmrCoreData.H"
#include "initialization/EmbeddedBoundary.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/spacecharge/ForceFromSelfFields.H"
#include "particles/spacecharge/GatherAndPush.H"
#include "particles/spacecharge/Gauss3dPush.H"
#include "particles/spacecharge/Gauss2p5dPush.H"
#include "particles/spacecharge/PoissonSolve.H"
#include "particles/transformation/CoordinateTransformation.H"

#include <AMReX_REAL.H>

#include <memory>
#include <stdexcept>


namespace impactx::particles::spacecharge
{
    void HandleSpacecharge (
        std::unique_ptr<initialization::AmrCoreData> & amr_data,
        std::function<void()> ResizeMesh,
        amrex::Real slice_ds
    )
    {
        if (slice_ds == 0.0) { return; }

        BL_PROFILE("impactx::particles::wakefields::HandleSpacecharge")

        auto space_charge = get_space_charge_algo();

        // turn off if disabled by user
        if (space_charge == SpaceChargeAlgo::False) { return; }

        // turn off if less than 2 particles
        if (amr_data->track_particles.m_particle_container->TotalNumberOfParticles(true, false) < 2) { return; }

        if (space_charge != SpaceChargeAlgo::True_2D && space_charge != SpaceChargeAlgo::True_2p5D && space_charge != SpaceChargeAlgo::Gauss2p5D)
        {
            // transform from x',y',t to x,y,z
            transformation::CoordinateTransformation(
                *amr_data->track_particles.m_particle_container,
                CoordSystem::t
            );
        }

        if (space_charge == SpaceChargeAlgo::Gauss3D)
        {
            Gauss3dPush(*amr_data->track_particles.m_particle_container, slice_ds);
        }
        else if (space_charge == SpaceChargeAlgo::Gauss2p5D)
        {
            Gauss2p5dPush(*amr_data->track_particles.m_particle_container, slice_ds);
        }
        else if (space_charge == SpaceChargeAlgo::True_3D || space_charge == SpaceChargeAlgo::True_2D || space_charge == SpaceChargeAlgo::True_2p5D)
        {
            // Note: The following operations assume that
            // the particles are in x, y, z coordinates.

            // Resize the mesh, based on `amr_data->track_particles.m_particle_container` extent
            ResizeMesh();

            // Redistribute particles in the new mesh in:
            // 3D: x, y, z
            // 2D: x, y, t
            amr_data->track_particles.m_particle_container->Redistribute();

            // charge deposition
            amr_data->track_particles.m_particle_container->DepositCharge(
                amr_data->track_particles.m_rho,
                amr_data->refRatio()
            );

            // poisson solve in x,y,z
            spacecharge::PoissonSolve(
                *amr_data->track_particles.m_particle_container,
                amr_data->track_particles.m_rho,
                amr_data->track_particles.m_rho_2d,
                amr_data->track_particles.m_phi,
                amr_data->refRatio(),
                initialization::eb_factories(*amr_data)
            );

            // calculate force in x,y,z
            spacecharge::ForceFromSelfFields(
                amr_data->track_particles.m_space_charge_field,
                amr_data->track_particles.m_phi,
                amr_data->Geom()
            );

            // gather and space-charge push in x,y,z , assuming the space-charge
            // field is the same before/after transformation
            // TODO: This is currently using linear order.
            spacecharge::GatherAndPush(
                *amr_data->track_particles.m_particle_container,
                amr_data->track_particles.m_space_charge_field,
                amr_data->track_particles.m_phi,
                amr_data->Geom(),
                slice_ds
            );
        }

        if (space_charge != SpaceChargeAlgo::True_2D && space_charge != SpaceChargeAlgo::True_2p5D && space_charge != SpaceChargeAlgo::Gauss2p5D)
        {
            // transform from x,y,z to x',y',t
            transformation::CoordinateTransformation(
                *amr_data->track_particles.m_particle_container,
                CoordSystem::s
            );
        }

        // for later: original Impact implementation as an option for 3D space charge to:
        // Redistribute particles in x',y',t
        //   TODO: only needed if we want to gather and push space charge
        //         in x',y',t
        //   TODO: change geometry beforehand according to transformation
        //m_amr_data->track_particles.m_particle_container->Redistribute();
        //
        // in original Impact, we gather and space-charge push in x',y',t ,
        // assuming that the distribution did not change
    }

} // namespace impactx::particles::spacecharge
