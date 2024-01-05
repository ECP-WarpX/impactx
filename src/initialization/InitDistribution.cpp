/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "initialization/InitDistribution.H"

#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/All.H"

#include <ablastr/constant.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <variant>


namespace impactx
{
    void
    ImpactX::add_particles (
        amrex::ParticleReal bunch_charge,
        distribution::KnownDistributions distr,
        int npart
    )
    {
        BL_PROFILE("ImpactX::add_particles");

        auto const & ref = m_particle_container->GetRefParticle();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ref.charge_qe() != 0.0,
            "add_particles: Reference particle charge not yet set!");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ref.mass_MeV() != 0.0,
            "add_particles: Reference particle mass not yet set!");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ref.kin_energy_MeV() != 0.0,
            "add_particles: Reference particle energy not yet set!");

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bunch_charge >= 0.0,
            "add_particles: the bunch charge should be positive. "
            "For negatively charge bunches, please change the reference particle's charge.");
        if (bunch_charge == 0.0) {
            ablastr::warn_manager::WMRecordWarning(
                "ImpactX::add_particles",
                "The bunch charge is set to zero. ImpactX will run with "
                "zero-weighted particles. Did you mean to set the space "
                "charge algorithm to off instead?",
                ablastr::warn_manager::WarnPriority::low
            );
        }

        // Logic: We initialize 1/Nth of particles, independent of their
        // position, per MPI rank. We then measure the distribution's spatial
        // extent, create a grid, resize it to fit the beam, and then
        // redistribute particles so that they reside on the correct MPI rank.
        int const myproc = amrex::ParallelDescriptor::MyProc();
        int const nprocs = amrex::ParallelDescriptor::NProcs();
        int const navg = npart / nprocs;
        int const nleft = npart - navg * nprocs;
        int npart_this_proc = (myproc < nleft) ? navg+1 : navg;
        auto const rel_part_this_proc =
            amrex::ParticleReal(npart_this_proc) / amrex::ParticleReal(npart);

        // alloc data for particle attributes
        amrex::Gpu::DeviceVector<amrex::ParticleReal> x, y, t;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> px, py, pt;
        x.resize(npart_this_proc);
        y.resize(npart_this_proc);
        t.resize(npart_this_proc);
        px.resize(npart_this_proc);
        py.resize(npart_this_proc);
        pt.resize(npart_this_proc);

        std::visit([&](auto&& distribution){
            // initialize distributions
            distribution.initialize(bunch_charge, ref);

            amrex::ParticleReal * const AMREX_RESTRICT x_ptr = x.data();
            amrex::ParticleReal * const AMREX_RESTRICT y_ptr = y.data();
            amrex::ParticleReal * const AMREX_RESTRICT t_ptr = t.data();
            amrex::ParticleReal * const AMREX_RESTRICT px_ptr = px.data();
            amrex::ParticleReal * const AMREX_RESTRICT py_ptr = py.data();
            amrex::ParticleReal * const AMREX_RESTRICT pt_ptr = pt.data();

            initialization::InitSingleParticleData const init_single_particle_data(
                distribution, x_ptr, y_ptr, t_ptr, px_ptr, py_ptr, pt_ptr);
            amrex::ParallelForRNG(npart_this_proc, init_single_particle_data);

            // finalize distributions and deallocate temporary device global memory
            amrex::Gpu::streamSynchronize();
            distribution.finalize();
        }, distr);

        m_particle_container->AddNParticles(x, y, t, px, py, pt,
                                            ref.qm_qeeV(),
                                            bunch_charge * rel_part_this_proc);

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
        amrex::ParmParse const pp_dist("beam");

        amrex::ParticleReal kin_energy = 0.0;  // Beam kinetic energy (MeV)
        pp_dist.get("kin_energy", kin_energy);

        amrex::ParticleReal bunch_charge = 0.0;  // Bunch charge (C)
        pp_dist.get("charge", bunch_charge);

        std::string particle_type;  // Particle type
        pp_dist.get("particle", particle_type);

        amrex::ParticleReal qe;     // charge (elementary charge)
        amrex::ParticleReal massE;  // MeV/c^2
        if (particle_type == "electron") {
            qe = -1.0;
            massE = ablastr::constant::SI::m_e / ablastr::constant::SI::MeV_invc2;
        } else if (particle_type == "positron") {
            qe = 1.0;
            massE = ablastr::constant::SI::m_e / ablastr::constant::SI::MeV_invc2;
        } else if (particle_type == "proton") {
            qe = 1.0;
            massE = ablastr::constant::SI::m_p / ablastr::constant::SI::MeV_invc2;
        }
        else {  // default to electron
            ablastr::warn_manager::WMRecordWarning(
                "ImpactX::initBeamDistributionFromInputs",
                "No beam.particle specified, defaulting to electrons.",
                ablastr::warn_manager::WarnPriority::low
            );
            qe = -1.0;
            massE = ablastr::constant::SI::m_e / ablastr::constant::SI::MeV_invc2;
        }

        // set charge and mass and energy of ref particle
        m_particle_container->GetRefParticle()
                .set_charge_qe(qe).set_mass_MeV(massE).set_kin_energy_MeV(kin_energy);

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

          distribution::KnownDistributions const waterbag(distribution::Waterbag(
              sigx, sigy, sigt,
              sigpx, sigpy, sigpt,
              muxpx, muypy, mutpt));

          add_particles(bunch_charge, waterbag, npart);

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

          distribution::KnownDistributions const kurth6D(distribution::Kurth6D(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(bunch_charge, kurth6D, npart);

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

          distribution::KnownDistributions const gaussian(distribution::Gaussian(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(bunch_charge, gaussian, npart);

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

          distribution::KnownDistributions const kvDist(distribution::KVdist(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(bunch_charge, kvDist, npart);

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

          distribution::KnownDistributions const kurth4D(distribution::Kurth4D(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(bunch_charge, kurth4D, npart);
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

          distribution::KnownDistributions const semigaussian(distribution::Semigaussian(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(bunch_charge, semigaussian, npart);

        } else if (distribution_type == "triangle") {
          amrex::ParticleReal sigx, sigy, sigt, sigpx, sigpy, sigpt;
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

          distribution::KnownDistributions const triangle(distribution::Triangle(
            sigx, sigy, sigt,
            sigpx, sigpy, sigpt,
            muxpx, muypy, mutpt));

          add_particles(bunch_charge, triangle, npart);

        } else if (distribution_type == "thermal") {
          amrex::ParticleReal k, kT, kT_halo, normalize, normalize_halo;
          amrex::ParticleReal halo = 0.0;
          pp_dist.get("k", k);
          pp_dist.get("kT", kT);
          kT_halo = kT;
          pp_dist.get("normalize", normalize);
          normalize_halo = normalize;
          pp_dist.query("kT_halo", kT_halo);
          pp_dist.query("normalize_halo", normalize_halo);
          pp_dist.query("halo", halo);

          distribution::KnownDistributions thermal(distribution::Thermal(k, kT, kT_halo, normalize, normalize_halo, halo));

          add_particles(bunch_charge, thermal, npart);

        } else {
            amrex::Abort("Unknown distribution: " + distribution_type);
        }

        // print information on the initialized beam
        amrex::Print() << "Beam kinetic energy (MeV): " << kin_energy << std::endl;
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
