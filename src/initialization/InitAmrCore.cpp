/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "InitAmrCore.H"

#include "initialization/InitAMReX.H"

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_CoordSys.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <stdexcept>


namespace impactx::initialization
{
namespace details
{

    /** Set the initial physical simulation extent
     *
     * Sets to either a fake min:max box or to the user-provided extent.
     * Controlled by geometry.dynamic_size.
     *
     * @return the RealBox in meters
     */
    amrex::RealBox
    init_physical_domain ()
    {
        amrex::ParmParse pp_geometry("geometry");
        bool dynamic_size = true;
        pp_geometry.queryAdd("dynamic_size", dynamic_size);
        amrex::Vector<amrex::Real> prob_lo(undefined_geometry_prob_lo.begin(), undefined_geometry_prob_lo.end());
        amrex::Vector<amrex::Real> prob_hi(undefined_geometry_prob_hi.begin(), undefined_geometry_prob_hi.end());
        if (dynamic_size)
        {
            pp_geometry.queryAdd("prob_lo", prob_lo);
            pp_geometry.queryAdd("prob_hi", prob_hi);
        }
        else
        {
            pp_geometry.getarr("prob_lo", prob_lo);
            pp_geometry.getarr("prob_hi", prob_hi);
        }
        amrex::RealBox rb(prob_lo.data(), prob_hi.data());
        return rb;
    }
} // namespace details

    AmrCoreData
    init_amr_core ()
    {
        if (!amrex::Initialized())
        {
            default_init_AMReX();
            // note: due to global state, it would be too early to use amrex::Abort() here
            //throw std::runtime_error("AMReX must be initialized before ImpactX simulation can be constructed.");
        }

        amrex::ParmParse pp_amr("amr");

        amrex::Vector<int> n_cell(AMREX_SPACEDIM);
        bool const has_ncell = pp_amr.queryarr("n_cell", n_cell);
        if (has_ncell)
            return amrex_amrcore_gridding();
        else
            return one_box_per_rank();
    }

    AmrCoreData
    amrex_amrcore_gridding ()
    {
        amrex::ParmParse pp_amr("amr");

        // Domain index space
        amrex::Vector<int> n_cell;
        pp_amr.getarr("n_cell", n_cell);
        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell));

        // Domain physical size
        //   we might resize this dynamically
        auto rb = details::init_physical_domain();

        // Periodicity (none)
        amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        // Mesh-refinement
        int max_level = 0;
        pp_amr.query("max_level", max_level);
        //   amrex::AmrMesh::InitAmrMesh will query amr.ref_ratio or amr.ref_ratio_vect on its own
        amrex::Vector<amrex::IntVect> const & ref_ratios = amrex::Vector<amrex::IntVect>();

        return {rb, max_level, n_cell, amrex::CoordSys::cartesian, ref_ratios, is_periodic};
    }

    AmrCoreData
    one_box_per_rank ()
    {
        // Domain index space
        amrex::AmrInfo amr_info;
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const amrex::IntVect high_end = amr_info.blocking_factor[0]
                                        * amrex::IntVect(AMREX_D_DECL(nprocs,1,1)) - amrex::IntVect(1);
        amrex::Box domain(amrex::IntVect(0), high_end);
        //   adding amr.n_cell for consistency
        amrex::ParmParse pp_amr("amr");
        auto const n_cell_iv = domain.size();
        amrex::Vector<int> n_cell_v(n_cell_iv.begin(), n_cell_iv.end());
        pp_amr.addarr("n_cell", n_cell_v);

        // Domain physical size
        //   we might resize this dynamically
        auto rb = details::init_physical_domain();

        // Periodicity (none)
        amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        amrex::Geometry geom(domain, rb, amrex::CoordSys::cartesian, is_periodic);
        return {geom, amr_info};
    }
} // namespace impactx::initialization
