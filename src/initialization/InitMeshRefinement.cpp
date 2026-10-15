/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/Algorithms.H"
#include "initialization/EmbeddedBoundary.H"
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/Waterbag.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx
{
namespace detail
{
    amrex::Vector<amrex::Real>
    read_mr_prob_relative ()
    {
        amrex::ParmParse pp_algo("algo");
        amrex::ParmParse pp_amr("amr");
        amrex::ParmParse pp_geometry("geometry");

        int max_level = 0;
        pp_amr.queryWithParser("max_level", max_level);

        std::string poisson_solver = "fft";
        pp_algo.queryAdd("poisson_solver", poisson_solver);

        // The box is expanded beyond the min and max of the particle beam.
        amrex::Vector<amrex::Real> prob_relative(max_level + 1, 1.0);
        prob_relative[0] = 3.0;  // top/bottom pad the beam on the lowest level by default by its width
        pp_geometry.queryarr("prob_relative", prob_relative);

        if (prob_relative[0] < 3.0 && poisson_solver == "multigrid")
            ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::read_mr_prob_relative",
                    "Dynamic resizing of the mesh uses a geometry.prob_relative "
                    "padding of less than 3 for level 0. This might result in boundary "
                    "artifacts for space charge calculation. "
                    "There is no minimum good value for this parameter, consider "
                    "doing a convergence test.",
                    ablastr::warn_manager::WarnPriority::high
            );

        if (prob_relative[0] < 1.0)
            throw std::runtime_error("geometry.prob_relative must be >= 1.0 (the beam size) on the coarsest level");

        // check that prob_relative[0] > prob_relative[1] > prob_relative[2] ...
        amrex::Real last_lev_rel = std::numeric_limits<amrex::Real>::max();
        for (int lev = 0; lev <= max_level; ++lev) {
            amrex::Real const prob_relative_lvl = prob_relative[lev];
            if (prob_relative_lvl <= 0.0)
                throw std::runtime_error("geometry.prob_relative must be strictly positive for all levels");
            if (prob_relative_lvl > last_lev_rel)
                throw std::runtime_error("geometry.prob_relative must be descending over refinement levels");

            last_lev_rel = prob_relative_lvl;
        }

        return prob_relative;
    }
}

    void ImpactX::ResizeMesh ()
    {
        BL_PROFILE("ImpactX::ResizeMesh");

        using namespace amrex::literals; // for _rt and _prt

        {
            auto space_charge = get_space_charge_algo();
            if (space_charge == SpaceChargeAlgo::False)
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::ResizeMesh",
                    "This is a simulation without space charge. "
                    "ResizeMesh (and pc.Redistribute) should only be called "
                    "in space charge simulations.",
                    ablastr::warn_manager::WarnPriority::high
                );
        }

        // Extract the min and max of the particle positions
        auto const [x_min, y_min, z_min, x_max, y_max, z_max] = amr_data->track_particles.m_particle_container->MinAndMaxPositions();

        // guard for flat beams:
        //   https://github.com/BLAST-ImpactX/impactx/issues/44
        if (x_min == x_max || y_min == y_max || z_min == z_max)
            throw std::runtime_error("Flat beam detected. This is not yet supported: https://github.com/BLAST-ImpactX/impactx/issues/44");

        auto const wall = initialization::read_transverse_wall();
        if (wall && (x_min < wall->bbox_lo[0] || x_max > wall->bbox_hi[0] ||
                     y_min < wall->bbox_lo[1] || y_max > wall->bbox_hi[1]))
        {
            ablastr::warn_manager::WMRecordWarning(
                "ImpactX::ResizeMesh",
                "Particles found outside of the eb.shape bounding box. "
                "Use element apertures to remove them.",
                ablastr::warn_manager::WarnPriority::medium
            );
        }

        amrex::ParmParse pp_geometry("geometry");
        bool dynamic_size = true;
        pp_geometry.query("dynamic_size", dynamic_size);

        amrex::Vector<amrex::RealBox> rb(amr_data->finestLevel() + 1);  // extent per level
        if (dynamic_size)
        {
            // The coarsest level is expanded (or reduced) relative the min and max of particles.
            auto const prob_relative = detail::read_mr_prob_relative();

            amrex::Real const frac = prob_relative[0];
            amrex::RealVect const beam_min(x_min, y_min, z_min);
            amrex::RealVect const beam_max(x_max, y_max, z_max);
            amrex::RealVect const beam_width(beam_max - beam_min);

            amrex::RealVect const beam_padding = beam_width * (frac - 1_rt) * 0.5_rt;
            //                           added to the beam extent --^         ^-- box half above/below the beam

            amrex::RealVect lo = beam_min - beam_padding;
            amrex::RealVect hi = beam_max + beam_padding;

            // transverse walls: x and y cover the aperture plus 2 cells on each side
            if (wall)
            {
                for (int dir = 0; dir < 2; ++dir)
                {
                    int const n = amr_data->Geom(0).Domain().length(dir);
                    if (n <= 4) {
                        throw std::runtime_error("eb.shape requires more than 4 cells in x and y");
                    }
                    amrex::Real const center = 0.5_rt * (wall->bbox_lo[dir] + wall->bbox_hi[dir]);
                    amrex::Real const width = wall->bbox_hi[dir] - wall->bbox_lo[dir];
                    amrex::Real const half_width =
                        0.5_rt * width * amrex::Real(n) / amrex::Real(n - 4);
                    lo[dir] = center - half_width;
                    hi[dir] = center + half_width;
                }
            }

            // In AMReX, all levels have the same problem domain, that of the
            // coarsest level, even if only partly covered.
            for (int lev = 0; lev <= amr_data->finestLevel(); ++lev)
            {
                rb[lev].setLo(lo);
                rb[lev].setHi(hi);
            }
        }
        else
        {
            // note: we read and set the size again because an interactive /
            //       Python user might have changed it between steps
            amrex::Vector<amrex::Real> prob_lo;
            amrex::Vector<amrex::Real> prob_hi;
            pp_geometry.getarr("prob_lo", prob_lo);
            pp_geometry.getarr("prob_hi", prob_hi);

            rb[0] = {prob_lo.data(), prob_hi.data()};

            if (wall && (prob_lo[0] > wall->bbox_lo[0] || prob_hi[0] < wall->bbox_hi[0] ||
                         prob_lo[1] > wall->bbox_lo[1] || prob_hi[1] < wall->bbox_hi[1]))
            {
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::ResizeMesh",
                    "geometry.prob_lo/hi does not contain the eb.shape aperture: "
                    "the domain faces cut through the beam pipe.",
                    ablastr::warn_manager::WarnPriority::high
                );
            }

            if (amr_data->maxLevel() > 1)
                amrex::Abort("Did not implement ResizeMesh for static domains and >1 MR levels.");
        }

        // updating geometry.prob_lo/hi for consistency
        amrex::Vector<amrex::Real> const prob_lo = {rb[0].lo()[0], rb[0].lo()[1], rb[0].lo()[2]};
        amrex::Vector<amrex::Real> const prob_hi = {rb[0].hi()[0], rb[0].hi()[1], rb[0].hi()[2]};
        pp_geometry.addarr("prob_lo", prob_lo);
        pp_geometry.addarr("prob_hi", prob_hi);

        // Resize the domain size
        amrex::Geometry::ResetDefaultProbDomain(rb[0]);

        for (int lev = 0; lev <= amr_data->finestLevel(); ++lev)
        {
            amrex::Geometry g = amr_data->Geom(lev);
            g.ProbDomain(rb[lev]);
            amr_data->SetGeometry(lev, g);

            amr_data->track_particles.m_particle_container->SetParticleGeometry(lev, g);
        }

        // the EB data depend on the geometry
        if (wall) {
            initialization::build_eb(*amr_data);
        }
    }
} // namespace impactx
