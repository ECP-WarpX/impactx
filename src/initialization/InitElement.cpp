/* Copyright 2021-2022 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/elements/All.H"

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <vector>


namespace impactx
{
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
} // namespace impactx
