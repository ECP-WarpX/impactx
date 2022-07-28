/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/All.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>


namespace impactx
{
    void
    ImpactX::add_particles (
        amrex::ParticleReal qm,
        amrex::ParticleReal bunch_charge,
        distribution::KnownDistributions distr,
        int npart
    )
    {
        BL_PROFILE("ImpactX::add_particles");

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

        std::visit([&](auto&& distribution){
            x.reserve(npart_this_proc);
            y.reserve(npart_this_proc);
            t.reserve(npart_this_proc);
            px.reserve(npart_this_proc);
            py.reserve(npart_this_proc);
            pt.reserve(npart_this_proc);

            for (amrex::Long i = 0; i < npart_this_proc; ++i) {
                distribution(ix, iy, it, ipx, ipy, ipt, rng);
                x.push_back(ix);
                y.push_back(iy);
                t.push_back(it);
                px.push_back(ipx);
                py.push_back(ipy);
                pt.push_back(ipt);
            }
        }, distr);

        int const lev = 0;
        m_particle_container->AddNParticles(lev, x, y, t, px, py, pt,
                                            qm, bunch_charge);

        // Resize the mesh to fit the spatial extent of the beam and then
        // redistribute particles, so they reside on the MPI rank that is
        // responsible for the respective spatial particle position.
        this->ResizeMesh();
        m_particle_container->Redistribute();
    }

    void ImpactX::initBeamDistributionFromInputs ()
    {
        BL_PROFILE("ImpactX::initBeamDistributionFromInputs");

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

          distribution::KnownDistributions waterbag(distribution::Waterbag(
              sigx, sigy, sigt,
              sigpx, sigpy, sigpt,
              muxpx, muypy, mutpt));

          add_particles(qm, bunch_charge, waterbag, npart);

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

          distribution::KnownDistributions kurth6D(distribution::Kurth6D(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(qm, bunch_charge, kurth6D, npart);

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

          distribution::KnownDistributions gaussian(distribution::Gaussian(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(qm, bunch_charge, gaussian, npart);

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

          distribution::KnownDistributions kvDist(distribution::KVdist(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(qm, bunch_charge, kvDist, npart);

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

          distribution::KnownDistributions kurth4D(distribution::Kurth4D(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(qm, bunch_charge, kurth4D, npart);
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

          distribution::KnownDistributions semigaussian(distribution::Semigaussian(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(qm, bunch_charge, semigaussian, npart);
        } else {
            amrex::Abort("Unknown distribution: " + distribution_type);
        }

        // reference particle
        amrex::ParticleReal massE;  // MeV
        if (particle_type == "electron") {
            massE = 0.510998950;
        } else if (particle_type == "proton") {
            massE = 938.27208816;
        } else {
            massE = 0.510998950;  // default to electron
        }

        // set the energy and initialize the reference particle
        RefPart refpart;
        refpart.Set_Energy_MeV(energy, massE);
        m_particle_container->SetRefParticle(refpart);

//        m_particle_container->Set_Energy_MeV(energy, massE);

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
        amrex::Print() << "# of particles: " << m_particle_container->TotalNumberOfParticles() << std::endl;
    }
} // namespace impactx
