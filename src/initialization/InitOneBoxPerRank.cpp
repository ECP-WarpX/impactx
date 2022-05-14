/* Copyright 2021-2022 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "InitOneBoxPerRank.H"

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_CoordSys.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>
#include <AMReX_SPACE.H>


namespace impactx::initialization
{
    AmrCoreData
    one_box_per_rank ()
    {
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
