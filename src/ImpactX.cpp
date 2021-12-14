/* Copyright 2021 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

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
        amrex::Print() << "boxArray(0) " << boxArray(0) << std::endl;;

        m_particle_container->AddNParticles(0, {0.0}, {0.2}, {0.4});
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
        // todo data_mf.define(ba, dm, 1, 0);
        amrex::ignore_unused(lev, time, ba, dm);
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
        // todo
        amrex::ignore_unused(lev);
    }

    void ImpactX::ResizeMesh () {
        // Extract the min and max of the particle positions
        auto const [x_min, x_max, y_min, y_max, z_min, z_max] = m_particle_container->MinAndMaxPositions();
        // Resize the domain size
        // The box is expanded slightly beyond the min and max of particles.
        // This controlled by the variable `frac` below.
        const amrex::Real frac=0.1;
        amrex::RealBox rb(
            {x_min-frac*(x_max-x_min), y_min-frac*(y_max-y_min), z_min-frac*(z_max-z_min)}, // Low bound
            {x_max+frac*(x_max-x_min), y_max+frac*(y_max-y_min), z_max+frac*(z_max-z_min)}); // High bound
        amrex::Geometry::ResetDefaultProbDomain(rb);
        for (int lev = 0; lev <= max_level; ++lev) {
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

            // Note: The following operation assume that
            // the particles are in x, y, z coordinates.
            // Resize the mesh, based on `m_particle_container` extent
            ResizeMesh();
            // Redistribute particles in the new mesh
            m_particle_container->Redistribute();

            // push all particles
            Push(*m_particle_container, m_lattice);

            // do more stuff in the step
            //...

            amrex::Print() << "\n";

        } // end step loop
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
        for (std::string const element_name : lattice_elements) {
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
            } else {
                amrex::Abort("Unknown type for lattice element " + element_name + ": " + element_type);
            }
        }

        amrex::Print() << "Initialized element list" << std::endl;
    }

} // namespace impactx
