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
           amrex::ParticleReal const k = 27.245985285371864_prt;
           amrex::ParticleReal const phi = -0.2_prt;
           amrex::ParticleReal const E0 = 62.0_prt;

           refpart.t = t;
           refpart.pt = pt - tau*E0*RF_Efield(zeval)*cos(k*t+phi);
           zeval = zeval;
           amrex::AllPrintToFile("rfout") << zeval << " " << RF_Efield(zeval) << "\n";
        }
    }

    amrex::ParticleReal
    RF_Efield (amrex::ParticleReal & zeval)
    {
        using namespace amrex::literals; // for _rt and _prt

        // specify constants
        amrex::ParticleReal const pi = 3.141592653589793_prt;
        amrex::ParticleReal const zlen = 1.31879807_prt;
        amrex::ParticleReal const zmid = zlen/2.0_prt;
        int const ncoef = 25;

        // specify Fourier coefficients
        amrex::ParticleReal cos_coef [ncoef] = {};
        amrex::ParticleReal sin_coef [ncoef] = {};

        // recall array indexing begins with 0
        cos_coef[0] = 0.1644024074311037;
        cos_coef[1] = -0.1324009958969339;
        cos_coef[2] = 4.3443060026047219E-002;
        cos_coef[3] = 8.5602654094946495E-002;
        cos_coef[4] = -0.2433578169042885;
        cos_coef[5] = 0.5297150596779437;
        cos_coef[6] = 0.7164884680963959;
        cos_coef[7] = -5.2579522442877296E-003;
        cos_coef[8] = -5.5025369142193678E-002;
        cos_coef[9] = 4.6845673335028933E-002;
        cos_coef[10] = -2.3279346335638568E-002;
        cos_coef[11] = 4.0800777539657775E-003;
        cos_coef[12] = 4.1378326533752169E-003;
        cos_coef[13] = -2.5040533340490805E-003;
        cos_coef[14] = -4.0654981400000964E-003;
        cos_coef[15] = 9.6630592067498289E-003;
        cos_coef[16] = -8.5275895985990214E-003;
        cos_coef[17] = -5.8078747006425020E-002;
        cos_coef[18] = -2.4044337836660403E-002;
        cos_coef[19] = 1.0968240064697212E-002;
        cos_coef[20] = -3.4461179858301418E-003;
        cos_coef[21] = -8.1201564869443749E-004;
        cos_coef[22] = 2.1438992904959380E-003;
        cos_coef[23] = -1.4997753525697276E-003;
        cos_coef[24] = 1.8685171825676386E-004;

        // compute on-axis electric field (z is relative to cavity midpoint)
        amrex::ParticleReal efield = 0.0;
        amrex::ParticleReal efieldp = 0.0;
        amrex::ParticleReal efieldpp = 0.0;
        amrex::ParticleReal const z = zeval - zmid;

        if (abs(z)<=zmid)
        {
           efield = 0.5_prt*cos_coef[0];
           for(int j=1; j < ncoef; ++j)
           {
             efield = efield + cos_coef[j]*cos(j*2*pi*z/zlen) + 
                  sin_coef[j]*sin(j*2*pi*z/zlen);
             efieldp = efieldp-j*2*pi*cos_coef[j]*sin(j*2*pi*z/zlen)/zlen +
                  j*2*pi*sin_coef[j]*cos(j*2*pi*z/zlen)/zlen;
             efieldpp = efieldpp- pow(j*2*pi*cos_coef[j]/zlen,2) *cos(j*2*pi*z/zlen) -
                  pow(j*2*pi*sin_coef[j]/zlen,2) *sin(j*2*pi*z/zlen);
           }
        }
        return efield;
    }
    
} // namespace impactx::integrators
