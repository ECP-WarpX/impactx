/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/Waterbag.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <string>
#include <vector>


namespace impactx
{
    /** Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
    {
        // todo
        amrex::ignore_unused(lev, tags, time, ngrow);
    }

    /** Make a new level from scratch using provided BoxArray and DistributionMapping.
     *
     * Only used during initialization.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm)
    {
        amrex::ignore_unused(time, ba, dm);

        // set human-readable tag for each MultiFab
        auto const tag = [lev]( std::string tagname ) {
            tagname.append("[l=").append(std::to_string(lev)).append("]");
            return amrex::MFInfo().SetTag(std::move(tagname));
        };

        // charge (rho) mesh
        amrex::BoxArray cba = ba;
        // for MR levels (TODO):
        //cba.coarsen(refRatio(lev - 1));

        // staggering and number of charge components in the field
        auto const rho_nodal_flag = amrex::IntVect::TheNodeVector();
        int const num_components_rho = 1;

        // guard cells for charge deposition
        int const particle_shape = m_particle_container->GetParticleShape();
        int num_guards_rho = 0;
        if (particle_shape % 2 == 0)  // even shape orders
            num_guards_rho = particle_shape / 2 + 1;
        else  // odd shape orders
            num_guards_rho = (particle_shape + 1) / 2;

        m_rho.emplace(
            lev,
            amrex::MultiFab{amrex::convert(cba, rho_nodal_flag), dm, num_components_rho, num_guards_rho, tag("rho")});

        // scalar potential
        auto const phi_nodal_flag = rho_nodal_flag;
        int const num_components_phi = 1;
        int const num_guards_phi = num_guards_rho + 1; // todo: I think this just depends on max(MLMG, force calc)
        m_phi.emplace(
            lev,
            amrex::MultiFab{amrex::convert(cba, phi_nodal_flag), dm, num_components_phi, num_guards_phi, tag("phi")});

        // space charge force
        for (std::string const comp : {"x", "y", "z"})
        {
            std::string const str_tag = "space_charge_force_" + comp;

            std::unordered_map<std::string, amrex::MultiFab> f_comp;
            f_comp.emplace(
                comp,
                amrex::MultiFab{
                    amrex::convert(cba, rho_nodal_flag),
                    dm,
                    num_components_rho,
                    num_guards_rho,
                    tag(str_tag)
                }
            );
            m_space_charge_force.emplace(lev, std::move(f_comp));
        }
    }

    /** Make a new level using provided BoxArray and DistributionMapping and fill
     *  with interpolated coarse level data.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm)
    {
        // todo
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /** Remake an existing level using provided BoxArray and DistributionMapping
     *  and fill with existing fine and coarse data.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm)
    {
        // todo
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /** Delete level data
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::ClearLevel (int lev)
    {
        m_rho.erase(lev);
        m_phi.erase(lev);
        m_space_charge_force.erase(lev);
    }

    void ImpactX::ResizeMesh ()
    {
        BL_PROFILE("ImpactX::ResizeMesh");

        // Extract the min and max of the particle positions
        auto const [x_min, y_min, z_min, x_max, y_max, z_max] = m_particle_container->MinAndMaxPositions();
        // Resize the domain size
        // The box is expanded slightly beyond the min and max of particles.
        // This controlled by the variable `frac` below.
        const amrex::Real frac=0.1;
        amrex::RealBox rb(
            {x_min-frac*(x_max-x_min), y_min-frac*(y_max-y_min), z_min-frac*(z_max-z_min)}, // Low bound
            {x_max+frac*(x_max-x_min), y_max+frac*(y_max-y_min), z_max+frac*(z_max-z_min)}); // High bound
        amrex::Geometry::ResetDefaultProbDomain(rb);
        for (int lev = 0; lev <= this->max_level; ++lev) {
            amrex::Geometry g = Geom(lev);
            g.ProbDomain(rb);
            amrex::AmrMesh::SetGeometry(lev, g);
        }
    }
} // namespace impactx
