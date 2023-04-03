/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/elements/All.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <vector>


namespace impactx
{
namespace detail
{    /** Resizing version of amrex::ParmParse::queryAdd
     *
     * Work-around for https://github.com/AMReX-Codes/amrex/pull/3220
     *
     * @tparam T vector type
     * @param[inout] pp the parameter parser on the element to query
     * @param[in] name key name
     * @param[inout] ref vector with default values
     * @return indicates if key existed previously
     */
    template <typename T>
    int queryAddResize (amrex::ParmParse& pp, const char* name, std::vector<T>& ref) {
        std::vector<T> empty;
        int exist = pp.queryarr(name, empty);
        if (exist) {
            ref.resize(empty.size());
            pp.queryarr(name, ref);
        }
        if (!exist && !ref.empty()) {
            pp.addarr(name, ref);
        }
        return exist;
    }
} // namespace detail

    void ImpactX::initLatticeElementsFromInputs ()
    {
        BL_PROFILE("ImpactX::initLatticeElementsFromInputs");

        // make sure the element sequence is empty
        m_lattice.clear();

        // Parse the lattice elements
        amrex::ParmParse pp_lattice("lattice");
        std::vector<std::string> lattice_elements;
        pp_lattice.queryarr("elements", lattice_elements);

        // Default number of slices per element
        int nslice_default = 1;
        pp_lattice.query("nslice", nslice_default);

        // Default number of map integration steps per slice
        int mapsteps_default = 10;  // used only in RF cavity

        // Loop through lattice elements
        for (std::string const & element_name : lattice_elements) {
            // Check the element type
            amrex::ParmParse pp_element(element_name);
            std::string element_type;
            pp_element.get("type", element_type);

            // Initialize the corresponding element according to its type
            if (element_type == "quad") {
                amrex::Real ds, k;
                int nslice = nslice_default;
                pp_element.get("ds", ds);
                pp_element.get("k", k);
                pp_element.queryAdd("nslice", nslice);
                m_lattice.emplace_back( Quad(ds, k, nslice) );
            } else if (element_type == "drift") {
                amrex::Real ds;
                int nslice = nslice_default;
                pp_element.get("ds", ds);
                pp_element.queryAdd("nslice", nslice);
                m_lattice.emplace_back( Drift(ds, nslice) );
            } else if (element_type == "sbend") {
                amrex::Real ds, rc;
                int nslice = nslice_default;
                pp_element.get("ds", ds);
                pp_element.get("rc", rc);
                pp_element.queryAdd("nslice", nslice);
                m_lattice.emplace_back( Sbend(ds, rc, nslice) );
            } else if (element_type == "dipedge") {
                amrex::Real psi, rc, g, K2;
                pp_element.get("psi", psi);
                pp_element.get("rc", rc);
                pp_element.get("g", g);
                pp_element.get("K2", K2);
                m_lattice.emplace_back( DipEdge(psi, rc, g, K2) );
            } else if (element_type == "constf") {
                amrex::Real ds, kx, ky, kt;
                int nslice = nslice_default;
                pp_element.get("ds", ds);
                pp_element.get("kx", kx);
                pp_element.get("ky", ky);
                pp_element.get("kt", kt);
                pp_element.queryAdd("nslice", nslice);
                m_lattice.emplace_back( ConstF(ds, kx, ky, kt, nslice) );
            } else if (element_type == "shortrf") {
                amrex::Real V, k;
                pp_element.get("V", V);
                pp_element.get("k", k);
                m_lattice.emplace_back( ShortRF(V, k) );
            } else if (element_type == "multipole") {
                int m;
                amrex::Real k_normal, k_skew;
                pp_element.get("multipole", m);
                pp_element.get("k_normal", k_normal);
                pp_element.get("k_skew", k_skew);
                m_lattice.emplace_back( Multipole(m, k_normal, k_skew) );
            } else if (element_type == "nonlinear_lens") {
                amrex::Real knll, cnll;
                pp_element.get("knll", knll);
                pp_element.get("cnll", cnll);
                m_lattice.emplace_back( NonlinearLens(knll, cnll) );
            } else if (element_type == "rfcavity") {
                amrex::Real ds, escale, freq, phase;
                int nslice = nslice_default;
                int mapsteps = mapsteps_default;
                RF_field_data ez;
                std::vector<amrex::ParticleReal> cos_coef = ez.default_cos_coef;
                std::vector<amrex::ParticleReal> sin_coef = ez.default_sin_coef;
                pp_element.get("ds", ds);
                pp_element.get("escale", escale);
                pp_element.get("freq", freq);
                pp_element.get("phase", phase);
                pp_element.queryAdd("mapsteps", mapsteps);
                pp_element.queryAdd("nslice", nslice);
                detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
                detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);
                m_lattice.emplace_back( RFCavity(ds, escale, freq, phase, cos_coef, sin_coef, mapsteps, nslice) );
            } else if (element_type == "solenoid") {
                amrex::Real ds, ks;
                int nslice = nslice_default;
                pp_element.get("ds", ds);
                pp_element.get("ks", ks);
                pp_element.queryAdd("nslice", nslice);
                m_lattice.emplace_back( Sol(ds, ks, nslice) );
            } else if (element_type == "prot") {
                amrex::Real phi_in,phi_out;
                pp_element.get("phi_in", phi_in);
                pp_element.get("phi_out", phi_out);
                phi_in = phi_in*ablastr::constant::math::pi/180.0;
                phi_out = phi_out*ablastr::constant::math::pi/180.0;
                m_lattice.emplace_back( PRot(phi_in, phi_out) );
            } else if (element_type == "solenoid_softedge") {
                amrex::Real ds, bscale;
                int nslice = nslice_default;
                int mapsteps = mapsteps_default;
                Sol_field_data bz;
                std::vector<amrex::ParticleReal> cos_coef = bz.default_cos_coef;
                std::vector<amrex::ParticleReal> sin_coef = bz.default_sin_coef;
                pp_element.get("ds", ds);
                pp_element.get("bscale", bscale);
                pp_element.queryAdd("mapsteps", mapsteps);
                pp_element.queryAdd("nslice", nslice);
                detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
                detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);
                m_lattice.emplace_back( SoftSolenoid(ds, bscale, cos_coef, sin_coef, mapsteps, nslice) );
            } else if (element_type == "quadrupole_softedge") {
                amrex::Real ds, gscale;
                int nslice = nslice_default;
                int mapsteps = mapsteps_default;
                Quad_field_data gz;
                std::vector<amrex::ParticleReal> cos_coef = gz.default_cos_coef;
                std::vector<amrex::ParticleReal> sin_coef = gz.default_sin_coef;
                pp_element.get("ds", ds);
                pp_element.get("gscale", gscale);
                pp_element.queryAdd("mapsteps", mapsteps);
                pp_element.queryAdd("nslice", nslice);
                detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
                detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);
                m_lattice.emplace_back( SoftQuadrupole(ds, gscale, cos_coef, sin_coef, mapsteps, nslice) );
            } else if (element_type == "beam_monitor") {
                std::string openpmd_name = element_name;
                pp_element.queryAdd("name", openpmd_name);
                std::string openpmd_backend = "default";
                pp_element.queryAdd("backend", openpmd_backend);
                std::string  openpmd_encoding {"g"};
                pp_element.queryAdd("encoding", openpmd_encoding);
                m_lattice.emplace_back( diagnostics::BeamMonitor(openpmd_name, openpmd_backend, openpmd_encoding) );
            } else {
                amrex::Abort("Unknown type for lattice element " + element_name + ": " + element_type);
            }
        }

        amrex::Print() << "Initialized element list" << std::endl;
    }
} // namespace impactx
