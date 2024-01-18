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

#include <algorithm>
#include <map>
#include <string>
#include <utility>
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
        int const exist = pp.queryarr(name, empty);
        if (exist) {
            ref.resize(empty.size());
            pp.queryarr(name, ref);
        }
        if (!exist && !ref.empty()) {
            pp.addarr(name, ref);
        }
        return exist;
    }

    /** Read the Thick element parameters ds and nslice
     *
     * @param pp_element the element being read
     * @param nslice_default the default number of slices to use if not specified
     * @return total element length (ds) and number of slices through it (nslice)
     */
    std::pair<amrex::ParticleReal, int>
    query_ds (amrex::ParmParse& pp_element, int nslice_default)
    {
        amrex::ParticleReal ds;
        int nslice = nslice_default;
        pp_element.get("ds", ds);
        pp_element.queryAdd("nslice", nslice);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nslice > 0,
                                         pp_element.getPrefix() + ".nslice must be > 0.");

        return {ds, nslice};
    }

    /** Read the Alignment parameters dx, dy and rotation from inputs
     *
     * @param pp_element the element being read
     * @return key-value pairs for dx, dy and rotation_degree
     */
    std::map<std::string, amrex::ParticleReal>
    query_alignment (amrex::ParmParse& pp_element)
    {
        amrex::ParticleReal dx = 0;
        amrex::ParticleReal dy = 0;
        amrex::ParticleReal rotation_degree = 0;
        pp_element.query("dx", dx);
        pp_element.query("dy", dy);
        pp_element.query("rotation", rotation_degree);

        std::map<std::string, amrex::ParticleReal> values = {
                {"dx", dx},
                {"dy", dy},
                {"rotation_degree", rotation_degree}
        };

        return values;
    }
} // namespace detail

    /** Read a lattice element
     *
     * Read a lattice element from amrex::ParmParse, initialize it and append it to m_lattice.
     *
     * @param[in] element_name element name
     * @param[inout] m_lattice the accelerator lattice
     * @param[in] nslice_default
     * @param[in] mapsteps_default
     */
    void read_element (std::string const & element_name,
                       std::list<KnownElements> & m_lattice,
                       int nslice_default,
                       int mapsteps_default)
    {
        // Check the element type
        amrex::ParmParse pp_element(element_name);
        std::string element_type;
        pp_element.get("type", element_type);

        // Initialize the corresponding element according to its type
        if (element_type == "quad")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal k;
            pp_element.get("k", k);

            m_lattice.emplace_back( Quad(ds, k, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "drift")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            m_lattice.emplace_back( Drift(ds, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "sbend")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal rc;
            pp_element.get("rc", rc);

            m_lattice.emplace_back( Sbend(ds, rc, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "cfbend")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal rc, k;
            pp_element.get("rc", rc);
            pp_element.get("k", k);

            m_lattice.emplace_back( CFbend(ds, rc, k, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "dipedge")
        {
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal psi, rc, g, K2;
            pp_element.get("psi", psi);
            pp_element.get("rc", rc);
            pp_element.get("g", g);
            pp_element.get("K2", K2);

            m_lattice.emplace_back( DipEdge(psi, rc, g, K2, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "constf")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::Real kx, ky, kt;
            pp_element.get("kx", kx);
            pp_element.get("ky", ky);
            pp_element.get("kt", kt);

            m_lattice.emplace_back( ConstF(ds, kx, ky, kt, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "buncher")
        {
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal V, k;
            pp_element.get("V", V);
            pp_element.get("k", k);

            m_lattice.emplace_back( Buncher(V, k, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "shortrf")
        {
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal V, freq;
            amrex::ParticleReal phase = -90.0;
            pp_element.get("V", V);
            pp_element.get("freq", freq);
            pp_element.queryAdd("phase", phase);

            m_lattice.emplace_back( ShortRF(V, freq, phase, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "multipole")
        {
            auto a = detail::query_alignment(pp_element);

            int m;
            amrex::ParticleReal k_normal, k_skew;
            pp_element.get("multipole", m);
            pp_element.get("k_normal", k_normal);
            pp_element.get("k_skew", k_skew);

            m_lattice.emplace_back( Multipole(m, k_normal, k_skew, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "nonlinear_lens")
        {
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal knll, cnll;
            pp_element.get("knll", knll);
            pp_element.get("cnll", cnll);

            m_lattice.emplace_back( NonlinearLens(knll, cnll, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "rfcavity")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal escale, freq, phase;
            int mapsteps = mapsteps_default;
            RF_field_data const ez;
            std::vector<amrex::ParticleReal> cos_coef = ez.default_cos_coef;
            std::vector<amrex::ParticleReal> sin_coef = ez.default_sin_coef;
            pp_element.get("escale", escale);
            pp_element.get("freq", freq);
            pp_element.get("phase", phase);
            pp_element.queryAdd("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
            detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);

            m_lattice.emplace_back( RFCavity(ds, escale, freq, phase, cos_coef, sin_coef, a["dx"], a["dy"], a["rotation_degree"], mapsteps, nslice) );
        } else if (element_type == "solenoid")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal ks;
            pp_element.get("ks", ks);

            m_lattice.emplace_back( Sol(ds, ks, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "prot")
        {
            amrex::ParticleReal phi_in, phi_out;
            pp_element.get("phi_in", phi_in);
            pp_element.get("phi_out", phi_out);

            m_lattice.emplace_back( PRot(phi_in, phi_out) );
        } else if (element_type == "solenoid_softedge")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal bscale;
            int mapsteps = mapsteps_default;
            Sol_field_data const bz;
            std::vector<amrex::ParticleReal> cos_coef = bz.default_cos_coef;
            std::vector<amrex::ParticleReal> sin_coef = bz.default_sin_coef;
            pp_element.get("bscale", bscale);
            pp_element.queryAdd("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
            detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);

            m_lattice.emplace_back( SoftSolenoid(ds, bscale, cos_coef, sin_coef, a["dx"], a["dy"], a["rotation_degree"], mapsteps, nslice) );
        } else if (element_type == "quadrupole_softedge")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal gscale;
            int mapsteps = mapsteps_default;
            Quad_field_data const gz;
            std::vector<amrex::ParticleReal> cos_coef = gz.default_cos_coef;
            std::vector<amrex::ParticleReal> sin_coef = gz.default_sin_coef;
            pp_element.get("gscale", gscale);
            pp_element.queryAdd("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
            detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);

            m_lattice.emplace_back( SoftQuadrupole(ds, gscale, cos_coef, sin_coef, a["dx"], a["dy"], a["rotation_degree"], mapsteps, nslice) );
        } else if (element_type == "drift_chromatic")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            m_lattice.emplace_back( ChrDrift(ds, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "quad_chromatic")
        {
            auto a = detail::query_alignment(pp_element);
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);

            amrex::ParticleReal k;
            int units = 0;
            pp_element.get("k", k);
            pp_element.queryAdd("units", units);

            m_lattice.emplace_back( ChrQuad(ds, k, units, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "drift_exact")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            m_lattice.emplace_back( ExactDrift(ds, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "sbend_exact")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal phi;
            amrex::ParticleReal B = 0.0;
            pp_element.get("phi", phi);
            pp_element.queryAdd("B", B);

            m_lattice.emplace_back( ExactSbend(ds, phi, B, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "uniform_acc_chromatic")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal ez, bz;
            pp_element.get("ez", ez);
            pp_element.get("bz", bz);

            m_lattice.emplace_back( ChrAcc(ds, ez, bz, a["dx"], a["dy"], a["rotation_degree"], nslice) );
        } else if (element_type == "thin_dipole")
        {
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal theta, rc;
            pp_element.get("theta", theta);
            pp_element.get("rc", rc);

            m_lattice.emplace_back( ThinDipole(theta, rc, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "kicker")
        {
            auto a = detail::query_alignment(pp_element);

            amrex::ParticleReal xkick, ykick;
            std::string units_str = "dimensionless";
            pp_element.get("xkick", xkick);
            pp_element.get("ykick", ykick);
            pp_element.queryAdd("units", units_str);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(units_str == "dimensionless" || units_str == "T-m",
                                             element_name + ".units must be \"dimensionless\" or \"T-m\"");
            Kicker::UnitSystem const units = units_str == "dimensionless" ?
                Kicker::UnitSystem::dimensionless :
                Kicker::UnitSystem::Tm;

            m_lattice.emplace_back( Kicker(xkick, ykick, units, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "aperture")
        {
            auto a = detail::query_alignment(pp_element);

            amrex::Real xmax, ymax;
            std::string shape_str = "rectangular";
            pp_element.get("xmax", xmax);
            pp_element.get("ymax", ymax);
            pp_element.queryAdd("shape", shape_str);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(shape_str == "rectangular" || shape_str == "elliptical",
                                             element_name + ".shape must be \"rectangular\" or \"elliptical\"");
            Aperture::Shape shape = shape_str == "rectangular" ?
                                        Aperture::Shape::rectangular :
                                        Aperture::Shape::elliptical;

            m_lattice.emplace_back( Aperture(xmax, ymax, shape, a["dx"], a["dy"], a["rotation_degree"]) );
        } else if (element_type == "beam_monitor")
        {
            std::string openpmd_name = element_name;
            pp_element.queryAdd("name", openpmd_name);
            std::string openpmd_backend = "default";
            pp_element.queryAdd("backend", openpmd_backend);
            std::string openpmd_encoding{"g"};
            pp_element.queryAdd("encoding", openpmd_encoding);

            m_lattice.emplace_back(diagnostics::BeamMonitor(openpmd_name, openpmd_backend, openpmd_encoding));
        } else if (element_type == "line")
        {
            // Parse the lattice elements for the sub-lattice in the line
            amrex::ParmParse pp_sub_lattice(element_name);
            std::vector<std::string> sub_lattice_elements;
            pp_sub_lattice.queryarr("elements", sub_lattice_elements);
            bool reverse = false;
            pp_sub_lattice.queryAdd("reverse", reverse);
            int repeat = 1;
            pp_sub_lattice.queryAdd("repeat", repeat);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(repeat >= 1,
                                             element_name + ".repeat must be >= 1");

            if (reverse)
                std::reverse(sub_lattice_elements.begin(), sub_lattice_elements.end());

            for (int n=0; n<repeat; ++n) {
                for (std::string const &sub_element_name: sub_lattice_elements) {
                    read_element(sub_element_name, m_lattice, nslice_default, mapsteps_default);
                }
            }
        } else {
            amrex::Abort("Unknown type for lattice element " + element_name + ": " + element_type);
        }
    }

    void ImpactX::initLatticeElementsFromInputs ()
    {
        BL_PROFILE("ImpactX::initLatticeElementsFromInputs");

        // make sure the element sequence is empty
        m_lattice.clear();

        amrex::ParmParse pp_lattice("lattice");

        // periods through the lattice
        int periods = 1;
        pp_lattice.queryAdd("periods", periods);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(periods >= 1,
                                         "lattice.periods must be >= 1");

        // Parse the lattice elements
        std::vector<std::string> lattice_elements;
        pp_lattice.queryarr("elements", lattice_elements);

        // reverse the lattice order
        bool reverse = false;
        pp_lattice.queryAdd("reverse", reverse);
        if (reverse)
            std::reverse(lattice_elements.begin(), lattice_elements.end());

        // Default number of slices per element
        int nslice_default = 1;
        pp_lattice.query("nslice", nslice_default);

        // Default number of map integration steps per slice
        int const mapsteps_default = 10;  // used only in RF cavity

        // Loop through lattice elements
        for (std::string const & element_name : lattice_elements) {
            read_element(element_name, m_lattice, nslice_default, mapsteps_default);
        }

        amrex::Print() << "Initialized element list" << std::endl;
    }
} // namespace impactx
