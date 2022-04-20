/* Copyright 2021-2022 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/Waterbag.H"
#include "particles/distribution/Kurth6D.H"
#include "particles/distribution/Gaussian.H"
#include "particles/distribution/KVdist.H"
#include "particles/distribution/Kurth4D.H"
#include "particles/distribution/Semigaussian.H"

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>


namespace
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
}

namespace impactx
{
    void ImpactX::initDist ()
    {
        using namespace amrex::literals;

        // Parse the beam distribution parameters
        amrex::ParmParse pp_dist("beam");

        amrex::ParticleReal energy = 0.0;  // Beam kinetic energy (MeV)
        pp_dist.get("energy", energy);

        amrex::ParticleReal bunch_charge = 0.0;  // Bunch charge (C)
        pp_dist.get("charge", bunch_charge);

        std::string particle_type;  // Particle type
        pp_dist.get("particle", particle_type);

        amrex::ParticleReal qm = 0.0; // charge/mass ratio (q_e/eV)
        if(particle_type == "electron") {
            qm = -1.0/0.510998950e6;
        } else if(particle_type == "proton") {
            qm = 1.0/938.27208816e6;
        }
        else {
            qm = 0.0;
        }

        int npart = 1;  // Number of simulation particles
        pp_dist.get("npart", npart);

        std::string unit_type;  // System of units
        pp_dist.get("units", unit_type);

        std::string distribution_type;  // Beam distribution type
        pp_dist.get("distribution", distribution_type);

        if(distribution_type == "waterbag"){
          amrex::ParticleReal sigx,sigy,sigt,sigpx,sigpy,sigpt;
          amrex::ParticleReal muxpx = 0.0, muypy = 0.0, mutpt = 0.0;
          pp_dist.get("sigmaX", sigx);
          pp_dist.get("sigmaY", sigy);
          pp_dist.get("sigmaT", sigt);
          pp_dist.get("sigmaPx", sigpx);
          pp_dist.get("sigmaPy", sigpy);
          pp_dist.get("sigmaPt", sigpt);
          pp_dist.query("muxpx", muxpx);
          pp_dist.query("muypy", muypy);
          pp_dist.query("mutpt", mutpt);

          impactx::distribution::Waterbag waterbag(
              sigx, sigy, sigt,
              sigpx, sigpy, sigpt,
              muxpx, muypy, mutpt);

          generate_add_particles(*m_particle_container, qm, bunch_charge, waterbag, npart);

        } else if (distribution_type == "kurth6d") {
          amrex::ParticleReal sigx,sigy,sigt,sigpx,sigpy,sigpt;
          amrex::ParticleReal muxpx = 0.0, muypy = 0.0, mutpt = 0.0;
          pp_dist.get("sigmaX", sigx);
          pp_dist.get("sigmaY", sigy);
          pp_dist.get("sigmaT", sigt);
          pp_dist.get("sigmaPx", sigpx);
          pp_dist.get("sigmaPy", sigpy);
          pp_dist.get("sigmaPt", sigpt);
          pp_dist.query("muxpx", muxpx);
          pp_dist.query("muypy", muypy);
          pp_dist.query("mutpt", mutpt);

          impactx::distribution::Kurth6D kurth6D(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt);

          generate_add_particles(*m_particle_container, qm, bunch_charge, kurth6D, npart);

        } else if (distribution_type == "gaussian") {
          amrex::ParticleReal sigx,sigy,sigt,sigpx,sigpy,sigpt;
          amrex::ParticleReal muxpx = 0.0, muypy = 0.0, mutpt = 0.0;
          pp_dist.get("sigmaX", sigx);
          pp_dist.get("sigmaY", sigy);
          pp_dist.get("sigmaT", sigt);
          pp_dist.get("sigmaPx", sigpx);
          pp_dist.get("sigmaPy", sigpy);
          pp_dist.get("sigmaPt", sigpt);
          pp_dist.query("muxpx", muxpx);
          pp_dist.query("muypy", muypy);
          pp_dist.query("mutpt", mutpt);

          impactx::distribution::Gaussian gaussian(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt);

          generate_add_particles(*m_particle_container, qm, bunch_charge, gaussian, npart);

        } else if (distribution_type == "kvdist") {
          amrex::ParticleReal sigx,sigy,sigt,sigpx,sigpy,sigpt;
          amrex::ParticleReal muxpx = 0.0, muypy = 0.0, mutpt = 0.0;
          pp_dist.get("sigmaX", sigx);
          pp_dist.get("sigmaY", sigy);
          pp_dist.get("sigmaT", sigt);
          pp_dist.get("sigmaPx", sigpx);
          pp_dist.get("sigmaPy", sigpy);
          pp_dist.get("sigmaPt", sigpt);
          pp_dist.query("muxpx", muxpx);
          pp_dist.query("muypy", muypy);
          pp_dist.query("mutpt", mutpt);

          impactx::distribution::KVdist kvDist(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt);

          generate_add_particles(*m_particle_container, qm, bunch_charge, kvDist, npart);

        } else if (distribution_type == "kurth4d") {
          amrex::ParticleReal sigx,sigy,sigt,sigpx,sigpy,sigpt;
          amrex::ParticleReal muxpx = 0.0, muypy = 0.0, mutpt = 0.0;
          pp_dist.get("sigmaX", sigx);
          pp_dist.get("sigmaY", sigy);
          pp_dist.get("sigmaT", sigt);
          pp_dist.get("sigmaPx", sigpx);
          pp_dist.get("sigmaPy", sigpy);
          pp_dist.get("sigmaPt", sigpt);
          pp_dist.query("muxpx", muxpx);
          pp_dist.query("muypy", muypy);
          pp_dist.query("mutpt", mutpt);

          impactx::distribution::Kurth4D kurth4D(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt);

          generate_add_particles(*m_particle_container, qm, bunch_charge, kurth4D, npart);
        } else if (distribution_type == "semigaussian") {
          amrex::ParticleReal sigx,sigy,sigt,sigpx,sigpy,sigpt;
          amrex::ParticleReal muxpx = 0.0, muypy = 0.0, mutpt = 0.0;
          pp_dist.get("sigmaX", sigx);
          pp_dist.get("sigmaY", sigy);
          pp_dist.get("sigmaT", sigt);
          pp_dist.get("sigmaPx", sigpx);
          pp_dist.get("sigmaPy", sigpy);
          pp_dist.get("sigmaPt", sigpt);
          pp_dist.query("muxpx", muxpx);
          pp_dist.query("muypy", muypy);
          pp_dist.query("mutpt", mutpt);

          impactx::distribution::Semigaussian semigaussian(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt);

          generate_add_particles(*m_particle_container, qm, bunch_charge, semigaussian, npart);
        } else {
            amrex::Abort("Unknown distribution: " + distribution_type);
        }

        // Resize the mesh to fit the spatial extent of the beam and then
        // redistribute particles, so they reside on the MPI rank that is
        // responsible for the respective spatial particle position.
        this->ResizeMesh();
        m_particle_container->Redistribute();

        // reference particle
        amrex::ParticleReal massE;  // MeV
        if (particle_type == "electron") {
            massE = 0.510998950;
        } else if (particle_type == "proton") {
            massE = 938.27208816;
        } else {
            massE = 0.510998950;  // default to electron
        }
        RefPart refPart;
        refPart.x = 0.0;
        refPart.y = 0.0;
        refPart.t = 0.0;
        refPart.px = 0.0;
        refPart.py = 0.0;
        refPart.pt = -energy/massE - 1.0_prt;
        m_particle_container->SetRefParticle(refPart);

        // print information on the initialized beam
        amrex::Print() << "Beam kinetic energy (MeV): " << energy << std::endl;
        amrex::Print() << "Bunch charge (C): " << bunch_charge << std::endl;
        amrex::Print() << "Particle type: " << particle_type << std::endl;
        amrex::Print() << "Number of particles: " << npart << std::endl;
        amrex::Print() << "Beam distribution type: " << distribution_type << std::endl;

        if (unit_type == "static") {
            amrex::Print() << "Static units" << std::endl;
        } else if (unit_type == "dynamic") {
            amrex::Print() << "Dynamic units" << std::endl;
        } else {
            amrex::Abort("Unknown units (static/dynamic): " + unit_type);
        }

        amrex::Print() << "Initialized beam distribution parameters" << std::endl;
    }
} // namespace impactx
