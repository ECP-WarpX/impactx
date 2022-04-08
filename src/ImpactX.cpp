/* Copyright 2021-2022 Axel Huebl, Chad Mitchell, Ji Qiang, Remi Lehe
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "particles/transformation/CoordinateTransformation.H"
#include "particles/diagnostics/DiagnosticOutput.H"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>


namespace impactx
{
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

    void ImpactX::evolve (int num_steps)
    {
        BL_PROFILE("ImpactX::evolve");

        // print initial particle distribution to file
        diagnostics::DiagnosticOutput(*m_particle_container,
                                      diagnostics::OutputType::PrintParticles,
                                      "diags/initial_beam.txt");

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
                                      diagnostics::OutputType::PrintParticles,
                                      "diags/output_beam.txt");
    }
} // namespace impactx
