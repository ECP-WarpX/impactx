/* Copyright 2021 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "particles/transformation/CoordinateTransformation.H"
#include "particles/distribution/Waterbag.H"
#include "particles/distribution/Kurth6D.H"
#include "particles/diagnostics/DiagnosticOutput.H"

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <string>
#include <vector>

namespace impactx
{
    void
    overwrite_amrex_parser_defaults()
    {
        amrex::ParmParse pp_amrex("amrex");

        // https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters
        bool abort_on_out_of_gpu_memory = true; // AMReX' default: false
        pp_amrex.query("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);
        pp_amrex.add("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);
    }

    ImpactX::ImpactX (amrex::Geometry const& simulation_geometry, amrex::AmrInfo const& amr_info)
        : AmrCore(simulation_geometry, amr_info),
          m_particle_container(std::make_unique<ImpactXParticleContainer>(this))
    {
    }

    void ImpactX::initData ()
    {
        AmrCore::InitFromScratch(0.0);
        amrex::Print() << "boxArray(0) " << boxArray(0) << std::endl;

        // move old diagnostics out of the way
        amrex::UtilCreateCleanDirectory("diags", true);

        this->initDist();
        amrex::Print() << "# of particles: " << m_particle_container->TotalNumberOfParticles() << std::endl;
    }

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
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm)
    {
        amrex::ignore_unused(time, ba, dm);

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
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::ClearLevel (int lev)
    {
        m_rho.erase(lev);
    }

    void ImpactX::ResizeMesh () {
        // Extract the min and max of the particle positions
        auto const [x_min, y_min, z_min, x_max, y_max, z_max] = m_particle_container->MinAndMaxPositions();
        // Resize the domain size
        // The box is expanded slightly beyond the min and max of particles.
        // This controlled by the variable `frac` below.
        const amrex::Real frac=0.1;
        amrex::RealBox rb(
            {x_min-frac*(x_max-x_min), y_min-frac*(y_max-y_min), z_min-frac*(z_max-z_min)}, // Low bound
            {x_max+frac*(x_max-x_min), y_max+frac*(y_max-y_min), z_max+frac*(z_max-z_min)}); // High bound
        amrex::Geometry::ResetDefaultProbDomain(rb);
        for (int lev = 0; lev <= this->max_level; ++lev) {
            amrex::Geometry g = Geom(lev);
            g.ProbDomain(rb);
            amrex::AmrMesh::SetGeometry(lev, g);
        }
    }

    void ImpactX::evolve (int num_steps)
    {
        BL_PROFILE("ImpactX::evolve");

        for (int step = 0; step < num_steps; ++step)
        {
            BL_PROFILE("ImpactX::evolve::step");
            amrex::Print() << " ++++ Starting step=" << step << "\n";

            // transform from x',y',t to x,y,z
            transformation::CoordinateTransformation(*m_particle_container,
                                                     transformation::Direction::T2Z);

            // Space-charge calculation: turn off if there is only 1 particle
            if (m_particle_container->TotalNumberOfParticles(false,false) > 1) {

                // Note: The following operation assume that
                // the particles are in x, y, z coordinates.

                // Resize the mesh, based on `m_particle_container` extent
                ResizeMesh();

                // Redistribute particles in the new mesh in x, y, z
                //m_particle_container->Redistribute();  // extra overload/arguments?

                // charge deposition
                m_particle_container->DepositCharge(m_rho, this->refRatio());

                // poisson solve in x,y,z
                //   TODO

                // gather and space-charge push in x,y,z , assuming the space-charge
                // field is the same before/after transformation
                //   TODO
            }

            // transform from x,y,z to x',y',t
            transformation::CoordinateTransformation(*m_particle_container,
                                                     transformation::Direction::Z2T);

            // for later: original Impact implementation as an option
            // Redistribute particles in x',y',t
            //   TODO: only needed if we want to gather and push space charge
            //         in x',y',t
            //   TODO: change geometry beforehand according to transformation
            //m_particle_container->Redistribute();
            //
            // in original Impact, we gather and space-charge push in x',y',t ,
            // assuming that the distribution did not change

            // push all particles with external maps
            Push(*m_particle_container, m_lattice);

            // just prints an empty newline at the end of the step
            amrex::Print() << "\n";

        } // end step loop

        // print final particle distribution to file
        diagnostics::DiagnosticOutput(*m_particle_container,
                                      diagnostics::OutputType::PrintParticles);

    }

    void ImpactX::initElements ()
    {
        // make sure the element sequence is empty
        m_lattice.clear();

        // Parse the lattice elements
        amrex::ParmParse pp_lattice("lattice");
        std::vector<std::string> lattice_elements;
        pp_lattice.queryarr("elements", lattice_elements);

        // Loop through lattice elements
        for (std::string const & element_name : lattice_elements) {
            // Check the element type
            amrex::ParmParse pp_element(element_name);
            std::string element_type;
            pp_element.get("type", element_type);
            // Initialize the corresponding element according to its type
            if (element_type == "quad") {
                amrex::Real ds, k;
                pp_element.get("ds", ds);
                pp_element.get("k", k);
                m_lattice.emplace_back( Quad(ds, k) );
            } else if (element_type == "drift") {
                amrex::Real ds;
                pp_element.get("ds", ds);
                m_lattice.emplace_back( Drift(ds) );
            } else if (element_type == "sbend") {
                amrex::Real ds, rc;
                pp_element.get("ds", ds);
                pp_element.get("rc", rc);
                m_lattice.emplace_back( Sbend(ds, rc) );
            } else if (element_type == "dipedge") {
                amrex::Real psi, rc, g, K2;
                pp_element.get("psi", psi);
                pp_element.get("rc", rc);
                pp_element.get("g", g);
                pp_element.get("K2", K2);
                m_lattice.emplace_back( DipEdge(psi, rc, g, K2) );
            } else if (element_type == "constf") {
                amrex::Real ds, kx, ky, kt;
                pp_element.get("ds", ds);
                pp_element.get("kx", kx);
                pp_element.get("ky", ky);
                pp_element.get("kt", kt);
                m_lattice.emplace_back( ConstF(ds, kx, ky, kt) );
            } else {
                amrex::Abort("Unknown type for lattice element " + element_name + ": " + element_type);
            }
        }

        amrex::Print() << "Initialized element list" << std::endl;
    }

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

        amrex::ParticleReal qm = 0.0; //charge/mass ratio (q_e/eV)
        if(particle_type == "electron"){
          qm = -1.0/0.510998950e6;
        } else if(particle_type == "proton"){
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

          impactx::distribution::Waterbag waterbag(sigx,sigy,sigt,sigpx,
                                 sigpy,sigpt,muxpx,muypy,mutpt);

          amrex::Vector<amrex::ParticleReal> x, y, t;
          amrex::Vector<amrex::ParticleReal> px, py, pt;
          amrex::RandomEngine rng;
          amrex::ParticleReal ix, iy, it, ipx, ipy, ipt;

          if (amrex::ParallelDescriptor::IOProcessor()) {
              x.reserve(npart);
              y.reserve(npart);
              t.reserve(npart);
              px.reserve(npart);
              py.reserve(npart);
              pt.reserve(npart);

              // write file header
              amrex::PrintToFile("diags/initial_beam.txt") << "x y t px py pt\n";

              for(amrex::Long i = 0; i < npart; ++i) {

                  waterbag(ix, iy, it, ipx, ipy, ipt, rng);
                  x.push_back(ix);
                  y.push_back(iy);
                  t.push_back(it);
                  px.push_back(ipx);
                  py.push_back(ipy);
                  pt.push_back(ipt);
                  amrex::PrintToFile("diags/initial_beam.txt")
                      << ix << " " << iy << " " << it << " "
                      << ipx << " " << ipy << " " << ipt << "\n";
              }
          }

          int const lev = 0;
          m_particle_container->AddNParticles(lev, x, y, t, px, py, pt,
                                              qm, bunch_charge);

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

          impactx::distribution::Kurth6D Kurth6D(sigx,sigy,sigt,sigpx,
                                 sigpy,sigpt,muxpx,muypy,mutpt);

          amrex::Vector<amrex::ParticleReal> x, y, t;
          amrex::Vector<amrex::ParticleReal> px, py, pt;
          amrex::RandomEngine rng;
          amrex::ParticleReal ix, iy, it, ipx, ipy, ipt;

          if (amrex::ParallelDescriptor::IOProcessor()) {
              x.reserve(npart);
              y.reserve(npart);
              t.reserve(npart);
              px.reserve(npart);
              py.reserve(npart);
              pt.reserve(npart);

              // write file header
              amrex::PrintToFile("diags/initial_beam.txt") << "x y t px py pt\n";

              for(amrex::Long i = 0; i < npart; ++i) {

                  Kurth6D(ix, iy, it, ipx, ipy, ipt, rng);
                  x.push_back(ix);
                  y.push_back(iy);
                  t.push_back(it);
                  px.push_back(ipx);
                  py.push_back(ipy);
                  pt.push_back(ipt);
                  amrex::PrintToFile("diags/initial_beam.txt")
                      << ix << " " << iy << " " << it << " "
                      << ipx << " " << ipy << " " << ipt << "\n";
              }
          }

          int const lev = 0;
          m_particle_container->AddNParticles(lev, x, y, t, px, py, pt,
                                              qm, bunch_charge);

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
