/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "InitOneBoxPerRank.H"

#include "initialization/InitAMReX.H"

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_CoordSys.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>
#include <AMReX_SPACE.H>

#include <stdexcept>


namespace impactx::initialization
{
    AmrCoreData
    one_box_per_rank ()
    {
        if (!amrex::Initialized())
        {
            default_init_AMReX();
            // note: due to global state, it would be too early to use amrex::Abort() here
            //throw std::runtime_error("AMReX must be initialized before ImpactX simulation can be constructed.");
        }

        amrex::AmrInfo amr_info;
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const amrex::IntVect high_end = amr_info.blocking_factor[0]
                                        * amrex::IntVect(AMREX_D_DECL(nprocs,1,1)) - amrex::IntVect(1);
        amrex::Box domain(amrex::IntVect(0), high_end); // Domain index space
        amrex::RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});// Domain physical size
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
        amrex::Geometry geom(domain, rb, amrex::CoordSys::cartesian, is_periodic);

        return {geom, amr_info};
    }
} // namespace impactx::initialization
