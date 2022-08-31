/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "Integrators.H"

#include <ablastr/particles/IndexHandling.H>

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal
#include <AMReX_Print.H>      // for PrintToFile


namespace impactx::integrators
{
    void symp4_integrate (RefPart & refpart,
                           amrex::ParticleReal const zin,
                           amrex::ParticleReal const zout,
                           int const nsteps, std::string element_type)
    {

        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const alpha = 1.0_prt - pow(2.0_prt,1.0/3.0);
        amrex::ParticleReal const tau2 = dz/(1.0_prt + alpha);
        amrex::ParticleReal const tau1 = tau2/2.0_prt;
        amrex::ParticleReal const tau3 = alpha*tau1;
        amrex::ParticleReal const tau4 = (alpha - 1.0_prt)*tau2;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            map1(tau1,refpart,zeval,element_type);
            map2(tau2,refpart,zeval,element_type);
            map1(tau3,refpart,zeval,element_type);
            map2(tau4,refpart,zeval,element_type);
            map1(tau3,refpart,zeval,element_type);
            map2(tau2,refpart,zeval,element_type);
            map1(tau1,refpart,zeval,element_type);
        }

    }

    void map1 (amrex::ParticleReal const tau, RefPart & refpart,
                           amrex::ParticleReal & zeval, std::string element_type)
    {

        using namespace amrex::literals; // for _rt and _prt

        if(element_type == "rfcavity"){

           amrex::ParticleReal const t = refpart.t;
           amrex::ParticleReal const pt = refpart.pt;

           if(pt < -1.0_prt){
              refpart.t = t + tau/sqrt(1.0_prt - pow(pt,-2));
              refpart.pt = pt;
           }
           else {
             refpart.t = t;
             refpart.pt = pt;
           }
           zeval = zeval + tau;
        }

    }

    void map2 (amrex::ParticleReal const tau, RefPart & refpart,
                           amrex::ParticleReal & zeval, std::string element_type)
    {

        using namespace amrex::literals; // for _rt and _prt

        if(element_type == "rfcavity"){

           amrex::ParticleReal const t = refpart.t;
           amrex::ParticleReal const pt = refpart.pt;
           amrex::ParticleReal const E0 = 1.0_prt;  //E0 = qEz/mc^2

           refpart.t = t;
           refpart.pt = pt - tau*E0;
           zeval = zeval;
        }
    }

} // namespace impactx::integrators
