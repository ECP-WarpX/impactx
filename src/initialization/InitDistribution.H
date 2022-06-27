/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "particles/ImpactXParticleContainer.H"

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_Vector.H>

#include <string>


namespace impactx
{
    /** Generate and add n particles to a particle container
     *
     * @tparam T_Distribution distribution function to draw from (type)
     * @param pc a particle container to add particles to
     * @param qm charge/mass ratio (q_e/eV)
     * @param bunch_charge bunch charge (C)
     * @param distr distribution function to draw from (object)
     * @param npart number of particles to draw
     */
    template <typename T_Distribution>
    void
    generate_add_particles (
        impactx::ImpactXParticleContainer & pc,
        amrex::ParticleReal qm,
        amrex::ParticleReal bunch_charge,
        T_Distribution & distr,
        int npart
    )
    {
        amrex::Vector<amrex::ParticleReal> x, y, t;
        amrex::Vector<amrex::ParticleReal> px, py, pt;
        amrex::ParticleReal ix, iy, it, ipx, ipy, ipt;
        amrex::RandomEngine rng;

        // Logic: We initialize 1/Nth of particles, independent of their
        // position, per MPI rank. We then measure the distribution's spatial
        // extent, create a grid, resize it to fit the beam, and then
        // redistribute particles so that they reside on the correct MPI rank.
        int myproc = amrex::ParallelDescriptor::MyProc();
        int nprocs = amrex::ParallelDescriptor::NProcs();
        int navg = npart / nprocs;
        int nleft = npart - navg * nprocs;
        int npart_this_proc = (myproc < nleft) ? navg+1 : navg;

        x.reserve(npart_this_proc);
        y.reserve(npart_this_proc);
        t.reserve(npart_this_proc);
        px.reserve(npart_this_proc);
        py.reserve(npart_this_proc);
        pt.reserve(npart_this_proc);

        for (amrex::Long i = 0; i < npart_this_proc; ++i) {
            distr(ix, iy, it, ipx, ipy, ipt, rng);
            x.push_back(ix);
            y.push_back(iy);
            t.push_back(it);
            px.push_back(ipx);
            py.push_back(ipy);
            pt.push_back(ipt);
        }

        int const lev = 0;
        pc.AddNParticles(lev, x, y, t, px, py, pt,
                         qm, bunch_charge);
    }
} // namespace impactx
