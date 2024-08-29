/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang, Marco Garten
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

#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
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

        auto const & ref = amr_data->m_particle_container->GetRefParticle();
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

            using Distribution = std::remove_reference_t< std::remove_cv_t<decltype(distribution)> >;
            initialization::InitSingleParticleData<Distribution> const init_single_particle_data(
                distribution, x_ptr, y_ptr, t_ptr, px_ptr, py_ptr, pt_ptr);
            amrex::ParallelForRNG(npart_this_proc, init_single_particle_data);

            // finalize distributions and deallocate temporary device global memory
            amrex::Gpu::streamSynchronize();
            distribution.finalize();
        }, distr);

        amr_data->m_particle_container->AddNParticles(x, y, t, px, py, pt,
                                                      ref.qm_ratio_SI(),
                                            bunch_charge * rel_part_this_proc);

        bool space_charge = false;
        amrex::ParmParse pp_algo("algo");
        pp_algo.queryAdd("space_charge", space_charge);

        // For pure tracking simulations, we keep the particles split equally
        // on all MPI ranks, and ignore spatial "RealBox" extents of grids.
        if (space_charge) {
            // Resize the mesh to fit the spatial extent of the beam and then
            // redistribute particles, so they reside on the MPI rank that is
            // responsible for the respective spatial particle position.
            this->ResizeMesh();
            amr_data->m_particle_container->Redistribute();
        }
    }

    void initialization::set_distribution_parameters_from_twiss_inputs (
        amrex::ParmParse const & pp_dist,
        amrex::ParticleReal& sigx, amrex::ParticleReal& sigy, amrex::ParticleReal& sigt,
        amrex::ParticleReal& sigpx, amrex::ParticleReal& sigpy, amrex::ParticleReal& sigpt,
        amrex::ParticleReal& muxpx, amrex::ParticleReal& muypy, amrex::ParticleReal& mutpt
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // Values to be read from input
        amrex::ParticleReal betax, betay, betat, emittx, emitty, emittt;
        // If alpha is zero the bunch is in focus
        amrex::ParticleReal alphax = 0.0, alphay = 0.0, alphat = 0.0;

        // Reading the input Twiss parameters
        pp_dist.query("alphaX", alphax);
        pp_dist.query("alphaY", alphay);
        pp_dist.query("alphaT", alphat);
        pp_dist.get("betaX", betax);
        pp_dist.get("betaY", betay);
        pp_dist.get("betaT", betat);
        pp_dist.get("emittX", emittx);
        pp_dist.get("emittY", emitty);
        pp_dist.get("emittT", emittt);

        if (betax <= 0.0_prt || betay <= 0.0_prt || betat <= 0.0_prt) {
            throw std::runtime_error("Input Error: The beta function values need to be non-zero positive values in all dimensions.");
        }

        if (emittx <= 0.0_prt || emitty <= 0.0_prt || emittt <= 0.0_prt) {
            throw std::runtime_error("Input Error: Emittance values need to be non-zero positive values in all dimensions.");
        }

        std::array<amrex::ParticleReal, 3> const alphas = {alphax, alphay, alphat};
        std::array<amrex::ParticleReal, 3> const betas = {betax, betay, betat};
        std::array<amrex::ParticleReal, 3> const emittances = {emittx, emitty, emittt};

        // calculate Twiss / Courant-Snyder gammas
        amrex::Vector<amrex::ParticleReal> gammas;
        for (size_t i = 0; i < alphas.size(); i++)
            gammas.push_back((1.0 + std::pow(alphas.at(i), 2)) / betas.at(i));

        amrex::Vector<amrex::ParticleReal> lambdas_pos;
        amrex::Vector<amrex::ParticleReal> lambdas_mom;
        amrex::Vector<amrex::ParticleReal> correlations;

        // calculate intersections of phase space ellipse with coordinate axes and the correlation factors
        for (size_t k = 0; k < betas.size(); k++){
            lambdas_pos.push_back(std::sqrt(emittances.at(k)/gammas.at(k)));
            lambdas_mom.push_back(std::sqrt(emittances.at(k)/betas.at(k)));

            correlations.push_back(alphas.at(k) / std::sqrt(betas.at(k) * gammas.at(k)));
        }

        sigx = lambdas_pos.at(0);
        sigy = lambdas_pos.at(1);
        sigt = lambdas_pos.at(2);
        sigpx = lambdas_mom.at(0);
        sigpy = lambdas_mom.at(1);
        sigpt = lambdas_mom.at(2);
        muxpx = correlations.at(0);
        muypy = correlations.at(1);
        mutpt = correlations.at(2);
    }

    void initialization::set_distribution_parameters_from_phase_space_inputs (
        amrex::ParmParse const & pp_dist,
        amrex::ParticleReal& sigx, amrex::ParticleReal& sigy, amrex::ParticleReal& sigt,
        amrex::ParticleReal& sigpx, amrex::ParticleReal& sigpy, amrex::ParticleReal& sigpt,
        amrex::ParticleReal& muxpx, amrex::ParticleReal& muypy, amrex::ParticleReal& mutpt
    )
    {
        pp_dist.get("lambdaX", sigx);
        pp_dist.get("lambdaY", sigy);
        pp_dist.get("lambdaT", sigt);
        pp_dist.get("lambdaPx", sigpx);
        pp_dist.get("lambdaPy", sigpy);
        pp_dist.get("lambdaPt", sigpt);
        pp_dist.query("muxpx", muxpx);
        pp_dist.query("muypy", muypy);
        pp_dist.query("mutpt", mutpt);
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
        amr_data->m_particle_container->GetRefParticle()
                .set_charge_qe(qe).set_mass_MeV(massE).set_kin_energy_MeV(kin_energy);

        int npart = 1;  // Number of simulation particles
        pp_dist.get("npart", npart);

        std::string unit_type;  // System of units
        pp_dist.get("units", unit_type);

        // Beam distributions that all share the same input signature for parameters from a beam ellipse
        std::set<std::string> distribution_types_from_beam_ellipse = {
                "gaussian", "kurth4d", "kurth6d", "kvdist", "semigaussian", "triangle", "waterbag"
        };

        std::string distribution_type;  // Beam distribution type
        pp_dist.get("distribution", distribution_type);

        std::string base_dist_type = distribution_type;
        // Position of the underscore for splitting off the suffix in case the distribution name either ends in "_from_twiss" or "_from_cs"
        std::size_t str_pos_from_twiss = distribution_type.rfind("_from_twiss");
        std::size_t str_pos_from_cs = distribution_type.rfind("_from_cs");
        bool initialize_from_twiss = false;

        if (str_pos_from_twiss != std::string::npos) { // "_from_twiss" is found
            // Calculate suffix and base type, consider length of "_from_twiss" = 12
            base_dist_type = distribution_type.substr(0, str_pos_from_twiss);
            initialize_from_twiss = true;
        } else if (str_pos_from_cs != std::string::npos) { // "_from_cs" is found
            // Calculate suffix and base type, consider length of "_from_cs" = 8
            base_dist_type = distribution_type.substr(0, str_pos_from_cs);
            initialize_from_twiss = true;
        }

        /* After separating a potential suffix from its base type, check if the base distribution type is in the set of
         * distributions that all share the same input signature.
         */
        if (distribution_types_from_beam_ellipse.find(base_dist_type) != distribution_types_from_beam_ellipse.end())
        {
            amrex::ParticleReal sigx, sigy, sigt, sigpx, sigpy, sigpt;
            amrex::ParticleReal muxpx = 0.0, muypy = 0.0, mutpt = 0.0;

            if (initialize_from_twiss)
            {
                initialization::set_distribution_parameters_from_twiss_inputs(
                        pp_dist,
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt
                );
            } else {
                initialization::set_distribution_parameters_from_phase_space_inputs(
                        pp_dist,
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt
                );
            }

            if(base_dist_type == "waterbag"){
                distribution::KnownDistributions const waterbag(distribution::Waterbag(
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt));

                add_particles(bunch_charge, waterbag, npart);
            } else if (base_dist_type == "kurth6d") {
                distribution::KnownDistributions const kurth6D(distribution::Kurth6D(
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt));

                add_particles(bunch_charge, kurth6D, npart);
            } else if (base_dist_type == "gaussian") {
                distribution::KnownDistributions const gaussian(distribution::Gaussian(
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt));

                add_particles(bunch_charge, gaussian, npart);
            } else if (base_dist_type == "kvdist") {
                distribution::KnownDistributions const kvDist(distribution::KVdist(
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt));

                add_particles(bunch_charge, kvDist, npart);
            } else if (base_dist_type == "kurth4d") {
                distribution::KnownDistributions const kurth4D(distribution::Kurth4D(
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt));

                add_particles(bunch_charge, kurth4D, npart);
            } else if (base_dist_type == "semigaussian") {
                distribution::KnownDistributions const semigaussian(distribution::Semigaussian(
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt));

                add_particles(bunch_charge, semigaussian, npart);
            } else if (base_dist_type == "triangle") {
                distribution::KnownDistributions const triangle(distribution::Triangle(
                        sigx, sigy, sigt,
                        sigpx, sigpy, sigpt,
                        muxpx, muypy, mutpt));

                add_particles(bunch_charge, triangle, npart);
            } else {
                throw std::runtime_error("Unknown distribution: " + distribution_type);
            }

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
            throw std::runtime_error("Unknown distribution: " + distribution_type);
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
            throw std::runtime_error("Unknown units (static/dynamic): " + unit_type);
        }

        amrex::Print() << "Initialized beam distribution parameters" << std::endl;
        amrex::Print() << "# of particles: " << amr_data->m_particle_container->TotalNumberOfParticles() << std::endl;
    }
} // namespace impactx
