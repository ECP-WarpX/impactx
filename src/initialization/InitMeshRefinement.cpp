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
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/Waterbag.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <stdexcept>
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
     */
    void ImpactX::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm)
    {
        amrex::ignore_unused(time);

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

        m_rho.emplace(lev,
                      amrex::MultiFab{amrex::convert(cba, rho_nodal_flag), dm, num_components_rho, num_guards_rho, tag("rho")});
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
     */
    void ImpactX::ClearLevel (int lev)
    {
        m_rho.erase(lev);
    }

    void ImpactX::ResizeMesh ()
    {
        BL_PROFILE("ImpactX::ResizeMesh");

        // Extract the min and max of the particle positions
        auto const [x_min, y_min, z_min, x_max, y_max, z_max] = m_particle_container->MinAndMaxPositions();

        // guard for flat beams:
        //   https://github.com/ECP-WarpX/impactx/issues/44
        if (x_min == x_max || y_min == y_max || z_min == z_max)
            throw std::runtime_error("Flat beam detected. This is not yet supported: https://github.com/ECP-WarpX/impactx/issues/44");

        amrex::ParmParse pp_geometry("geometry");
        bool dynamic_size = true;
        pp_geometry.query("dynamic_size", dynamic_size);

        amrex::RealBox rb;
        if (dynamic_size)
        {
            // The box is expanded beyond the min and max of particles.
            // This controlled by the variable `frac` below.
            amrex::Real frac = 1.0;
            pp_geometry.query("prob_relative", frac);

            if (frac < 1.0)
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::ResizeMesh",
                    "Dynamic resizing of the mesh uses a geometry.prob_relative "
                    "padding of less than 1. This might result in boundary "
                    "artifacts for space charge calculation. "
                    "There is no minimum good value for this parameter, consider "
                    "doing a convergence test.",
                    ablastr::warn_manager::WarnPriority::high
                );

            if (frac < 0.0)
                throw std::runtime_error("geometry.prob_relative must be positive");

            rb = {
                {x_min - frac * (x_max - x_min), y_min - frac * (y_max - y_min),
                 z_min - frac * (z_max - z_min)}, // Low bound
                {x_max + frac * (x_max - x_min), y_max + frac * (y_max - y_min),
                 z_max + frac * (z_max - z_min)}  // High bound
            };
        }
        else
        {
            // note: we read and set the size again because an interactive /
            //       Python user might have changed it between steps
            amrex::Vector<amrex::Real> prob_lo;
            amrex::Vector<amrex::Real> prob_hi;
            pp_geometry.getarr("prob_lo", prob_lo);
            pp_geometry.getarr("prob_hi", prob_hi);

            rb = {prob_lo.data(), prob_hi.data()};
        }

        // updating geometry.prob_lo/hi for consistency
        amrex::Vector<amrex::Real> const prob_lo = {rb.lo()[0], rb.lo()[1], rb.lo()[2]};
        amrex::Vector<amrex::Real> const prob_hi = {rb.hi()[0], rb.hi()[1], rb.hi()[2]};
        pp_geometry.addarr("prob_lo", prob_lo);
        pp_geometry.addarr("prob_hi", prob_hi);

        // Resize the domain size
        amrex::Geometry::ResetDefaultProbDomain(rb);
        for (int lev = 0; lev <= this->max_level; ++lev) {
            amrex::Geometry g = Geom(lev);
            g.ProbDomain(rb);
            amrex::AmrMesh::SetGeometry(lev, g);
        }
    }
} // namespace impactx
