/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "elements/All.H"
#include "elements/mixin/lineartransport.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Enum.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>


namespace impactx
{
namespace detail
{
    /** Resizing version of amrex::ParmParse::queryAdd
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
        pp_element.getWithParser("ds", ds);
        pp_element.queryAddWithParser("nslice", nslice);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nslice > 0,
                                         pp_element.getPrefix() + ".nslice must be > 0.");

        return {ds, nslice};
    }

    /** Read the Alignment parameters dx, dy and rotation from inputs
     *
     * @param pp_element the element being read
     * @return key-value pairs for dx, dy and rotation_degree
     */
    template <typename T_Element>
    std::map<std::string, amrex::ParticleReal>
    query_alignment (amrex::ParmParse& pp_element)
    {
        amrex::ParticleReal dx = T_Element::DEFAULT_dx;
        amrex::ParticleReal dy = T_Element::DEFAULT_dy;
        amrex::ParticleReal rotation_degree = T_Element::DEFAULT_rotation_degree;
        pp_element.queryWithParser("dx", dx);
        pp_element.queryWithParser("dy", dy);
        pp_element.queryWithParser("rotation", rotation_degree);

        std::map<std::string, amrex::ParticleReal> values = {
                {"dx", dx},
                {"dy", dy},
                {"rotation_degree", rotation_degree}
        };

        return values;
    }

    /** Read the Aperture parameters aperture_x and aperture_y from inputs
     *
     * @param pp_element the element being read
     * @return key-value pairs for aperture_x and aperture_y
     */
    template <typename T_Element>
    std::map<std::string, amrex::ParticleReal>
    query_aperture (amrex::ParmParse& pp_element)
    {
        amrex::ParticleReal aperture_x = T_Element::DEFAULT_aperture_x;
        amrex::ParticleReal aperture_y = T_Element::DEFAULT_aperture_y;
        pp_element.query("aperture_x", aperture_x);
        pp_element.query("aperture_y", aperture_y);

        std::map<std::string, amrex::ParticleReal> values = {
                {"aperture_x", aperture_x},
                {"aperture_y", aperture_y}
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
     */
    void read_element (std::string const & element_name,
                       std::list<elements::KnownElements> & m_lattice,
                       int nslice_default)
    {
        using namespace amrex::literals; // for _prt
        using namespace elements;
;
        // Check the element type
        amrex::ParmParse pp_element(element_name);
        std::string element_type;
        pp_element.get("type", element_type);

        // Initialize the corresponding element according to its type
        if (element_type == "quad")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<Quad>(pp_element);
            auto b = detail::query_aperture<Quad>(pp_element);

            amrex::ParticleReal k;
            pp_element.getWithParser("k", k);

            m_lattice.emplace_back( Quad(ds, k, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "drift")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<Drift>(pp_element);
            auto b = detail::query_aperture<Drift>(pp_element);

            m_lattice.emplace_back( Drift(ds, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "sbend")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<Sbend>(pp_element);
            auto b = detail::query_aperture<Sbend>(pp_element);

            amrex::ParticleReal rc;
            pp_element.getWithParser("rc", rc);

            m_lattice.emplace_back( Sbend(ds, rc, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "cfbend")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<CFbend>(pp_element);
            auto b = detail::query_aperture<CFbend>(pp_element);

            amrex::ParticleReal rc, k;
            pp_element.getWithParser("rc", rc);
            pp_element.getWithParser("k", k);

            m_lattice.emplace_back( CFbend(ds, rc, k, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "dipedge")
        {
            auto a = detail::query_alignment<DipEdge>(pp_element);
            auto b = detail::query_aperture<DipEdge>(pp_element);

            amrex::ParticleReal psi, rc, g;
            amrex::ParticleReal R = DipEdge::DEFAULT_R;
            std::string model_str = amrex::getEnumNameString(DipEdge::DEFAULT_model);
            std::string location_str = amrex::getEnumNameString(DipEdge::DEFAULT_location);
            bool modify_ref_part = DipEdge::DEFAULT_modify_ref_part;

            amrex::ParticleReal K0 = DipEdge::DEFAULT_K0;
            amrex::ParticleReal K1 = DipEdge::DEFAULT_K1;
            amrex::ParticleReal K2 = DipEdge::DEFAULT_K2;
            amrex::ParticleReal K3 = DipEdge::DEFAULT_K3;
            amrex::ParticleReal K4 = DipEdge::DEFAULT_K4;
            amrex::ParticleReal K5 = DipEdge::DEFAULT_K5;
            amrex::ParticleReal K6 = DipEdge::DEFAULT_K6;
            pp_element.getWithParser("psi", psi);
            pp_element.getWithParser("rc", rc);
            pp_element.getWithParser("g", g);
            pp_element.queryAddWithParser("R", R);
            pp_element.queryAddWithParser("K0", K0);
            pp_element.queryAddWithParser("K1", K1);
            pp_element.queryAddWithParser("K2", K2);
            pp_element.queryAddWithParser("K3", K3);
            pp_element.queryAddWithParser("K4", K4);
            pp_element.queryAddWithParser("K5", K5);
            pp_element.queryAddWithParser("K6", K6);

            pp_element.queryAdd("model", model_str);
            DipEdge::Model const model = amrex::getEnum<DipEdge::Model>(model_str);
            pp_element.queryAdd("location", location_str);
            DipEdge::Location const location = amrex::getEnum<DipEdge::Location>(location_str);
            pp_element.queryAdd("modify_ref_part", modify_ref_part);

            if (R <= 0) {
                throw std::runtime_error(element_name + ".R must be >0 but is: " + std::to_string(R));
            }

            m_lattice.emplace_back( DipEdge(psi, rc, g, R, K0, K1, K2, K3, K4, K5, K6, model, location, modify_ref_part, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "quadedge")
        {
            auto a = detail::query_alignment<QuadEdge>(pp_element);
            auto b = detail::query_aperture<QuadEdge>(pp_element);

            amrex::ParticleReal k;
            int units = QuadEdge::DEFAULT_unit;
            std::string flag_str = amrex::getEnumNameString(QuadEdge::DEFAULT_flag);
            pp_element.getWithParser("k", k);
            pp_element.queryAddWithParser("units", units);
            pp_element.queryAdd("flag", flag_str);
            QuadEdge::Location const flag = amrex::getEnum<QuadEdge::Location>(flag_str);

            m_lattice.emplace_back( QuadEdge(k, units, flag, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );

        } else if (element_type == "constf")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ConstF>(pp_element);
            auto b = detail::query_aperture<ConstF>(pp_element);

            amrex::Real kx, ky, kt;
            pp_element.getWithParser("kx", kx);
            pp_element.getWithParser("ky", ky);
            pp_element.getWithParser("kt", kt);

            m_lattice.emplace_back( ConstF(ds, kx, ky, kt, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "buncher")
        {
            auto a = detail::query_alignment<Buncher>(pp_element);
            auto b = detail::query_aperture<Buncher>(pp_element);

            amrex::ParticleReal V, k;
            pp_element.getWithParser("V", V);
            pp_element.getWithParser("k", k);

            m_lattice.emplace_back( Buncher(V, k, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "shortrf")
        {
            auto a = detail::query_alignment<ShortRF>(pp_element);
            auto b = detail::query_aperture<ShortRF>(pp_element);

            amrex::ParticleReal V, freq;
            amrex::ParticleReal phase = ShortRF::DEFAULT_phase;
            pp_element.getWithParser("V", V);
            pp_element.getWithParser("freq", freq);
            pp_element.queryAddWithParser("phase", phase);

            m_lattice.emplace_back( ShortRF(V, freq, phase, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "multipole")
        {
            auto a = detail::query_alignment<Multipole>(pp_element);
            auto b = detail::query_aperture<Multipole>(pp_element);

            int m;
            amrex::ParticleReal k_normal, k_skew;
            pp_element.getWithParser("multipole", m);
            pp_element.getWithParser("k_normal", k_normal);
            pp_element.getWithParser("k_skew", k_skew);

            m_lattice.emplace_back( Multipole(m, k_normal, k_skew, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "nonlinear_lens")
        {
            auto a = detail::query_alignment<NonlinearLens>(pp_element);
            auto b = detail::query_aperture<NonlinearLens>(pp_element);

            amrex::ParticleReal knll, cnll;
            pp_element.getWithParser("knll", knll);
            pp_element.getWithParser("cnll", cnll);

            m_lattice.emplace_back( NonlinearLens(knll, cnll, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "rfcavity")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<RFCavity>(pp_element);
            auto b = detail::query_aperture<RFCavity>(pp_element);

            amrex::ParticleReal escale, freq, phase;
            int mapsteps = RFCavity::DEFAULT_mapsteps;
            RF_field_data const ez;
            std::vector<amrex::ParticleReal> cos_coef = ez.default_cos_coef;
            std::vector<amrex::ParticleReal> sin_coef = ez.default_sin_coef;
            std::vector<amrex::ParticleReal> z_data;
            std::vector<amrex::ParticleReal> wake_x_data;
            std::vector<amrex::ParticleReal> wake_y_data;
            std::vector<amrex::ParticleReal> wake_z_data;
            pp_element.getWithParser("escale", escale);
            pp_element.getWithParser("freq", freq);
            pp_element.getWithParser("phase", phase);
            pp_element.queryAddWithParser("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
            detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);
            detail::queryAddResize(pp_element, "z_data", z_data);
            detail::queryAddResize(pp_element, "wake_x_data", wake_x_data);
            detail::queryAddResize(pp_element, "wake_y_data", wake_y_data);
            detail::queryAddResize(pp_element, "wake_z_data", wake_z_data);

            m_lattice.emplace_back( RFCavity(ds, escale, freq, phase, cos_coef, sin_coef, z_data, wake_x_data, wake_y_data, wake_z_data, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], mapsteps, nslice, element_name) );
        } else if (element_type == "solenoid")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<Sol>(pp_element);
            auto b = detail::query_aperture<Sol>(pp_element);

            amrex::ParticleReal ks;
            pp_element.getWithParser("ks", ks);

            m_lattice.emplace_back( Sol(ds, ks, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "prot")
        {
            amrex::ParticleReal phi_in, phi_out;
            pp_element.getWithParser("phi_in", phi_in);
            pp_element.getWithParser("phi_out", phi_out);

            m_lattice.emplace_back( PRot(phi_in, phi_out, element_name) );
        } else if (element_type == "plane_xyrotation")
        {
            auto a = detail::query_alignment<PlaneXYRot>(pp_element);

            amrex::ParticleReal phi;
            pp_element.getWithParser("angle", phi);

            m_lattice.emplace_back( PlaneXYRot(phi, a["dx"], a["dy"], a["rotation_degree"], element_name) );
        } else if (element_type == "solenoid_softedge")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<SoftSolenoid>(pp_element);
            auto b = detail::query_aperture<SoftSolenoid>(pp_element);

            amrex::ParticleReal bscale;
            int mapsteps = SoftSolenoid::DEFAULT_mapsteps;
            int units = SoftSolenoid::DEFAULT_unit;
            Sol_field_data const bz;
            std::vector<amrex::ParticleReal> cos_coef = bz.default_cos_coef;
            std::vector<amrex::ParticleReal> sin_coef = bz.default_sin_coef;
            pp_element.getWithParser("bscale", bscale);
            pp_element.queryAddWithParser("units", units);
            pp_element.queryAddWithParser("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
            detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);

            m_lattice.emplace_back( SoftSolenoid(ds, bscale, cos_coef, sin_coef, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], mapsteps, nslice, element_name) );
        } else if (element_type == "quadrupole_softedge")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<SoftQuadrupole>(pp_element);
            auto b = detail::query_aperture<SoftQuadrupole>(pp_element);

            amrex::ParticleReal gscale;
            int mapsteps = SoftQuadrupole::DEFAULT_mapsteps;
            Quad_field_data const gz;
            std::vector<amrex::ParticleReal> cos_coef = gz.default_cos_coef;
            std::vector<amrex::ParticleReal> sin_coef = gz.default_sin_coef;
            pp_element.getWithParser("gscale", gscale);
            pp_element.queryAddWithParser("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "cos_coefficients", cos_coef);
            detail::queryAddResize(pp_element, "sin_coefficients", sin_coef);

            m_lattice.emplace_back( SoftQuadrupole(ds, gscale, cos_coef, sin_coef, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], mapsteps, nslice, element_name) );
        } else if (element_type == "drift_chromatic")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ChrDrift>(pp_element);
            auto b = detail::query_aperture<ChrDrift>(pp_element);

            m_lattice.emplace_back( ChrDrift(ds, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "quad_chromatic")
        {
            auto a = detail::query_alignment<ChrQuad>(pp_element);
            auto b = detail::query_aperture<ChrQuad>(pp_element);
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);

            amrex::ParticleReal k;
            int units = ChrQuad::DEFAULT_unit;
            pp_element.getWithParser("k", k);
            pp_element.queryAddWithParser("units", units);

            m_lattice.emplace_back( ChrQuad(ds, k, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "plasma_lens_chromatic")
        {
            auto a = detail::query_alignment<ChrPlasmaLens>(pp_element);
            auto b = detail::query_aperture<ChrPlasmaLens>(pp_element);
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);

            amrex::ParticleReal k;
            int units = ChrPlasmaLens::DEFAULT_unit;
            pp_element.getWithParser("k", k);
            pp_element.queryAddWithParser("units", units);

            m_lattice.emplace_back( ChrPlasmaLens(ds, k, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "tapered_plasma_lens")
        {
            auto a = detail::query_alignment<TaperedPL>(pp_element);
            auto b = detail::query_aperture<TaperedPL>(pp_element);

            amrex::ParticleReal k;
            amrex::ParticleReal taper;
            int units = TaperedPL::DEFAULT_unit;
            pp_element.getWithParser("k", k);
            pp_element.getWithParser("taper", taper);
            pp_element.queryAddWithParser("units", units);

            m_lattice.emplace_back( TaperedPL(k, taper, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "drift_exact")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ExactDrift>(pp_element);
            auto b = detail::query_aperture<ExactDrift>(pp_element);

            m_lattice.emplace_back( ExactDrift(ds, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "quad_exact")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ExactQuad>(pp_element);
            auto b = detail::query_aperture<ExactQuad>(pp_element);

            amrex::ParticleReal k;
            int units = ExactQuad::DEFAULT_unit;
            int int_order = ExactQuad::DEFAULT_int_order;
            int mapsteps = ExactQuad::DEFAULT_mapsteps;
            pp_element.getWithParser("k", k);
            pp_element.queryAddWithParser("units", units);
            pp_element.queryAddWithParser("int_order", int_order);
            pp_element.queryAddWithParser("mapsteps", mapsteps);

            m_lattice.emplace_back( ExactQuad(ds, k, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], int_order, mapsteps, nslice,
element_name) );
        } else if (element_type == "multipole_exact")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ExactMultipole>(pp_element);
            auto b = detail::query_aperture<ExactMultipole>(pp_element);

            std::vector<amrex::ParticleReal> k_normal = {0.0};
            std::vector<amrex::ParticleReal> k_skew = {0.0};
            int units = ExactMultipole::DEFAULT_unit;
            int int_order = ExactMultipole::DEFAULT_int_order;
            int mapsteps = ExactMultipole::DEFAULT_mapsteps;
            pp_element.queryAddWithParser("units", units);
            pp_element.queryAddWithParser("int_order", int_order);
            pp_element.queryAddWithParser("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "k_normal", k_normal);
            detail::queryAddResize(pp_element, "k_skew", k_skew);

            m_lattice.emplace_back( ExactMultipole(ds, k_normal, k_skew, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], int_order, mapsteps, nslice, element_name) );
        } else if (element_type == "cfbend_exact")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ExactCFbend>(pp_element);
            auto b = detail::query_aperture<ExactCFbend>(pp_element);

            std::vector<amrex::ParticleReal> k_normal = {0.0};
            std::vector<amrex::ParticleReal> k_skew = {0.0};
            int units = ExactCFbend::DEFAULT_unit;
            int int_order = ExactCFbend::DEFAULT_int_order;
            int mapsteps = ExactCFbend::DEFAULT_mapsteps;
            pp_element.queryAddWithParser("units", units);
            pp_element.queryAddWithParser("int_order", int_order);
            pp_element.queryAddWithParser("mapsteps", mapsteps);
            detail::queryAddResize(pp_element, "k_normal", k_normal);
            detail::queryAddResize(pp_element, "k_skew", k_skew);

            m_lattice.emplace_back( ExactCFbend(ds, k_normal, k_skew, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], int_order, mapsteps, nslice, element_name) );
        } else if (element_type == "sbend_exact")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ExactSbend>(pp_element);
            auto b = detail::query_aperture<ExactSbend>(pp_element);

            amrex::ParticleReal phi;
            amrex::ParticleReal B = ExactSbend::DEFAULT_B;
            pp_element.getWithParser("phi", phi);
            pp_element.queryAddWithParser("B", B);

            m_lattice.emplace_back( ExactSbend(ds, phi, B, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "uniform_acc_chromatic")
        {
            auto const [ds, nslice] = detail::query_ds(pp_element, nslice_default);
            auto a = detail::query_alignment<ChrAcc>(pp_element);
            auto b = detail::query_aperture<ChrAcc>(pp_element);

            amrex::ParticleReal ez, bz;
            pp_element.getWithParser("ez", ez);
            pp_element.getWithParser("bz", bz);

            m_lattice.emplace_back( ChrAcc(ds, ez, bz, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], nslice, element_name) );
        } else if (element_type == "thin_dipole")
        {
            auto a = detail::query_alignment<ThinDipole>(pp_element);
            auto b = detail::query_aperture<ThinDipole>(pp_element);

            amrex::ParticleReal theta, rc;
            pp_element.getWithParser("theta", theta);
            pp_element.getWithParser("rc", rc);

            m_lattice.emplace_back( ThinDipole(theta, rc, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "kicker")
        {
            auto a = detail::query_alignment<Kicker>(pp_element);
            auto b = detail::query_aperture<Kicker>(pp_element);

            amrex::ParticleReal xkick, ykick;
            std::string units_str = Kicker::unit_name(Kicker::DEFAULT_unit);
            pp_element.getWithParser("xkick", xkick);
            pp_element.getWithParser("ykick", ykick);
            pp_element.queryAdd("units", units_str);
            Kicker::UnitSystem const units = Kicker::unit_from_name(units_str);

            m_lattice.emplace_back( Kicker(xkick, ykick, units, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"], element_name) );
        } else if (element_type == "aperture")
        {
            auto a = detail::query_alignment<Aperture>(pp_element);

            amrex::Real aperture_x, aperture_y;
            amrex::ParticleReal repeat_x = Aperture::DEFAULT_repeat_x;
            amrex::ParticleReal repeat_y = Aperture::DEFAULT_repeat_y;
            bool shift_odd_x = Aperture::DEFAULT_shift_odd_x;
            std::string shape_str = amrex::getEnumNameString(Aperture::DEFAULT_shape);
            std::string action_str = amrex::getEnumNameString(Aperture::DEFAULT_action);

            // In the future, just use this:
            // pp_element.getWithParser("aperture_x", aperture_x);
            // pp_element.getWithParser("aperture_y", aperture_y);
            // Backwards compatibility to ImpactX <= 25.01
            bool const has_old_xmax = pp_element.queryWithParser("xmax", aperture_x);
            bool const has_old_ymax = pp_element.queryWithParser("ymax", aperture_y);
            if (has_old_xmax) {
                pp_element.queryAddWithParser("aperture_x", aperture_x);
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::read_element",
                    element_name + ".xmax is deprecated. Use " + element_name + ".aperture_x instead.",
                    ablastr::warn_manager::WarnPriority::high
                );
            } else {
                pp_element.getWithParser("aperture_x", aperture_x);
            }
            if (has_old_ymax) {
                pp_element.queryAddWithParser("aperture_y", aperture_y);
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::read_element",
                    element_name + ".ymax is deprecated. Use " + element_name + ".aperture_y instead.",
                    ablastr::warn_manager::WarnPriority::high
                );
            } else {
                pp_element.getWithParser("aperture_y", aperture_y);
            }

            pp_element.queryAddWithParser("repeat_x", repeat_x);
            pp_element.queryAddWithParser("repeat_y", repeat_y);
            pp_element.queryAdd("shift_odd_x", shift_odd_x);  // https://github.com/AMReX-Codes/amrex/issues/4535
            pp_element.queryAdd("shape", shape_str);
            pp_element.queryAdd("action", action_str);
            Aperture::Shape shape = amrex::getEnum<Aperture::Shape>(shape_str);
            Aperture::Action action = amrex::getEnum<Aperture::Action>(action_str);

            m_lattice.emplace_back( Aperture(aperture_x, aperture_y, repeat_x, repeat_y, shift_odd_x, shape, action, a["dx"], a["dy"], a["rotation_degree"], element_name) );
        } else if (element_type == "beam_monitor")
        {
            std::string openpmd_name = element_name;
            pp_element.queryAdd("name", openpmd_name);
            using impactx::elements::diagnostics::BeamMonitor;

            std::string openpmd_backend = BeamMonitor::DEFAULT_backend;
            pp_element.queryAdd("backend", openpmd_backend);
            std::string openpmd_encoding{BeamMonitor::DEFAULT_encoding};
            pp_element.queryAdd("encoding", openpmd_encoding);
            int period_sample_intervals = BeamMonitor::DEFAULT_period_sample_intervals;
            pp_element.queryAddWithParser("period_sample_intervals", period_sample_intervals);
            bool particles = BeamMonitor::DEFAULT_particles;
            pp_element.queryAdd("particles", particles);
            bool beam_moments = BeamMonitor::DEFAULT_beam_moments;
            pp_element.queryAdd("beam_moments", beam_moments);

            // optional: add and calculate additional particle properties
            // property: nonlinear lens invariants
            bool add_nll_invariants = false;
            pp_element.queryAdd("nonlinear_lens_invariants", add_nll_invariants);
            if (add_nll_invariants)
            {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(particles,
                    element_name + ".nonlinear_lens_invariants=true requires " + element_name + ".particles=true");

                amrex::ParticleReal alpha = 0.0;
                pp_element.queryAddWithParser("alpha", alpha);
                amrex::ParticleReal beta = 1.0;
                pp_element.queryAddWithParser("beta", beta);
                amrex::ParticleReal tn = 0.4_prt;
                pp_element.queryAddWithParser("tn", tn);
                amrex::ParticleReal cn = 0.01_prt;
                pp_element.queryAddWithParser("cn", cn);
            }

            m_lattice.emplace_back(BeamMonitor(openpmd_name, openpmd_backend, openpmd_encoding, period_sample_intervals, particles, beam_moments));
        } else if (element_type == "source")
        {
            std::string distribution, openpmd_path;
            pp_element.get("distribution", distribution);

            if (distribution == "openPMD")
            {
                pp_element.get("openpmd_path", openpmd_path);
            }

            bool active_once = Source::DEFAULT_active_once;
            pp_element.queryAdd("active_once", active_once);

            bool load_ref_particle = Source::DEFAULT_load_ref_particle;
            pp_element.queryAdd("load_ref_particle", load_ref_particle);

            // at most one of the two is set, the default reads the last step in the file
            std::optional<int> load_step = Source::DEFAULT_load_step;
            int load_step_value = 0;
            if (pp_element.queryWithParser("load_step", load_step_value)) {
                load_step = load_step_value;
            }

            std::optional<int> load_step_index = Source::DEFAULT_load_step_index;
            int load_step_index_value = 0;
            if (pp_element.queryWithParser("load_step_index", load_step_index_value)) {
                load_step_index = load_step_index_value;
            }

            m_lattice.emplace_back( Source(distribution, openpmd_path, active_once, load_ref_particle, load_step, load_step_index, element_name) );
        } else if (element_type == "line")
        {
            // Parse the lattice elements for the sub-lattice in the line
            amrex::ParmParse pp_sub_lattice(element_name);
            std::vector<std::string> sub_lattice_elements;
            pp_sub_lattice.queryarr("elements", sub_lattice_elements);
            bool reverse = false;
            pp_sub_lattice.queryAdd("reverse", reverse);
            int repeat = 1;
            pp_sub_lattice.queryAddWithParser("repeat", repeat);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(repeat >= 1,
                                             element_name + ".repeat must be >= 1");

            if (reverse)
                std::reverse(sub_lattice_elements.begin(), sub_lattice_elements.end());

            for (int n=0; n<repeat; ++n) {
                for (std::string const &sub_element_name: sub_lattice_elements) {
                    read_element(sub_element_name, m_lattice, nslice_default);
                }
            }
        } else if (element_type == "linear_map")
        {
            auto a = detail::query_alignment<LinearMap>(pp_element);
            auto b = detail::query_aperture<LinearMap>(pp_element);

            amrex::ParticleReal ds = LinearMap::DEFAULT_ds;
            pp_element.queryAdd("ds", ds);

            Map6x6 transport_map = Map6x6::Identity();

            // safe to ParmParse inputs for reproducibility
            for (int i=1; i<=6; ++i) {
                for (int j=1; j<=6; ++j) {
                    std::string name = "R" + std::to_string(i) + std::to_string(j);
                    pp_element.queryAddWithParser<amrex::ParticleReal>(name.c_str(), transport_map(i, j));
                }
            }

            m_lattice.emplace_back(LinearMap(transport_map, ds, a["dx"], a["dy"], a["rotation_degree"], b["aperture_x"], b["aperture_y"]) );
        } else if (element_type == "polygon_aperture")
        {
            auto a = detail::query_alignment<PolygonAperture>(pp_element);

            amrex::ParticleReal repeat_x = PolygonAperture::DEFAULT_repeat_x;
            amrex::ParticleReal repeat_y = PolygonAperture::DEFAULT_repeat_y;
            bool shift_odd_x = PolygonAperture::DEFAULT_shift_odd_x;
            std::string action_str = amrex::getEnumNameString(PolygonAperture::DEFAULT_action);

            std::vector<amrex::ParticleReal> vertices_x = {0.0};
            std::vector<amrex::ParticleReal> vertices_y = {0.0};
            amrex::ParticleReal min_radius2 = PolygonAperture::DEFAULT_min_radius2;

            detail::queryAddResize(pp_element, "vertices_x", vertices_x);
            detail::queryAddResize(pp_element, "vertices_y", vertices_y);
            pp_element.queryAddWithParser("min_radius2", min_radius2);

            pp_element.queryAddWithParser("repeat_x", repeat_x);
            pp_element.queryAddWithParser("repeat_y", repeat_y);
            pp_element.queryAdd("shift_odd_x", shift_odd_x);  // https://github.com/AMReX-Codes/amrex/issues/4535
            pp_element.queryAdd("action", action_str);

            //AMREX_ALWAYS_ASSERT_WITH_MESSAGE(action_str == "transmit" || action_str == "absorb",
            //                                 element_name + ".action must be \"transmit\" or \"absorb\"");

            PolygonAperture::Action action = amrex::getEnum<PolygonAperture::Action>(action_str);

            m_lattice.emplace_back(PolygonAperture(vertices_x, vertices_y, min_radius2, repeat_x, repeat_y,
                shift_odd_x, action, a["dx"], a["dy"], a["rotation_degree"], element_name) );

        } else if (element_type == "spin_map")
        {
            auto a = detail::query_alignment<SpinMap>(pp_element);

            amrex::ParticleReal ds = SpinMap::DEFAULT_ds;
            pp_element.queryAdd("ds", ds);

            Vector3 spin_rotation_vector = {};
            Map3x6 spin_orbit_coupling = {};

            // ParmParse inputs for spin rotation vector
            for (int i=1; i<=3; ++i) {
                 std::string name = "v" + std::to_string(i);
                 pp_element.queryAddWithParser<amrex::ParticleReal>(name.c_str(), spin_rotation_vector(i, 1));
            }

            // ParmParse inputs for spin-orbit coupling matrix
            for (int i=1; i<=3; ++i) {
                for (int j=1; j<=6; ++j) {
                    std::string name = "A" + std::to_string(i) + std::to_string(j);
                    pp_element.queryAddWithParser<amrex::ParticleReal>(name.c_str(), spin_orbit_coupling(i, j));
                }
            }

            m_lattice.emplace_back(SpinMap(spin_rotation_vector, spin_orbit_coupling, ds, a["dx"], a["dy"], a["rotation_degree"]) );

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
        pp_lattice.queryAddWithParser("periods", periods);
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
        int nslice_default = elements::Drift::DEFAULT_nslice;
        pp_lattice.queryWithParser("nslice", nslice_default);

        // Loop through lattice elements
        for (std::string const & element_name : lattice_elements) {
            read_element(element_name, m_lattice, nslice_default);
        }

        amrex::Print() << "Initialized element list" << std::endl;
    }
} // namespace impactx
