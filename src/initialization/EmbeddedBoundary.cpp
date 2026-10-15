/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Weiqun Zhang
 * License: BSD-3-Clause-LBNL
 */
#include "EmbeddedBoundary.H"

#include "initialization/AmrCoreData.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EB2.H>
#   include <AMReX_EB2_IF_Parser.H>
#   include <AMReX_EBFabFactory.H>
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>


namespace impactx::initialization
{
namespace
{
    /** Format a number for a parser expression, at full precision */
    std::string
    num (amrex::Real v)
    {
        std::ostringstream ss;
        ss << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10) << "(" << v << ")";
        return ss.str();
    }

    std::vector<amrex::Real>
    get_array (amrex::ParmParse const & pp, std::string const & name)
    {
        int const n = pp.countval(name);
        if (n == 0) {
            throw std::runtime_error("eb." + name + " is required for eb.shape");
        }
        std::vector<amrex::Real> v(n);
        pp.getarrWithParser(name, n, v.data());
        return v;
    }

    /** Max over the edges of a convex polygon of the outward half-plane distance */
    std::string
    convex_polygon_function (std::vector<amrex::Real> vx, std::vector<amrex::Real> vy)
    {
        auto const n = static_cast<int>(vx.size());
        if (n < 3 || vy.size() != vx.size()) {
            throw std::runtime_error(
                "eb.vertices_x and eb.vertices_y must have the same length of at least 3");
        }

        // orient counter-clockwise
        amrex::Real area2 = 0;
        for (int i = 0; i < n; ++i) {
            int const j = (i + 1) % n;
            area2 += vx[i] * vy[j] - vx[j] * vy[i];
        }
        if (area2 == 0) {
            throw std::runtime_error("eb.shape = polygon: the polygon has zero area");
        }
        if (area2 < 0) {
            std::reverse(vx.begin(), vx.end());
            std::reverse(vy.begin(), vy.end());
        }

        std::string f;
        for (int i = 0; i < n; ++i) {
            int const j = (i + 1) % n;
            int const k = (i + 2) % n;
            amrex::Real const dx = vx[j] - vx[i];
            amrex::Real const dy = vy[j] - vy[i];
            amrex::Real const len = std::hypot(dx, dy);
            if (len == 0) {
                throw std::runtime_error("eb.shape = polygon: repeated vertex");
            }
            // left turns only
            amrex::Real const cross = dx * (vy[k] - vy[j]) - dy * (vx[k] - vx[j]);
            if (cross <= 0) {
                throw std::runtime_error("eb.shape = polygon: the polygon must be convex");
            }
            amrex::Real const nx = dy / len;
            amrex::Real const ny = -dx / len;
            amrex::Real const c = nx * vx[i] + ny * vy[i];
            std::string const fi = num(nx) + "*x+" + num(ny) + "*y-" + num(c);
            f = f.empty() ? fi : "max(" + f + "," + fi + ")";
        }
        return f;
    }
} // namespace

    std::optional<TransverseWall>
    read_transverse_wall ()
    {
        amrex::ParmParse pp_eb("eb");
        std::string shape = "none";
        pp_eb.queryAdd("shape", shape);

        if (shape == "none") {
            return std::nullopt;
        }

        TransverseWall wall;
        if (shape == "elliptical" || shape == "rectangular")
        {
            amrex::Real a = 0, b = 0;
            pp_eb.getWithParser("aperture_x", a);
            pp_eb.getWithParser("aperture_y", b);
            if (a <= 0 || b <= 0) {
                throw std::runtime_error("eb.aperture_x and eb.aperture_y must be positive");
            }
            wall.implicit_function = (shape == "elliptical")
                ? "(x/" + num(a) + ")^2+(y/" + num(b) + ")^2-1"
                : "max(abs(x)/" + num(a) + ",abs(y)/" + num(b) + ")-1";
            wall.bbox_lo = {-a, -b};
            wall.bbox_hi = {a, b};
        }
        else if (shape == "polygon")
        {
            auto const vx = get_array(pp_eb, "vertices_x");
            auto const vy = get_array(pp_eb, "vertices_y");
            wall.implicit_function = convex_polygon_function(vx, vy);
            wall.bbox_lo = {*std::min_element(vx.begin(), vx.end()),
                            *std::min_element(vy.begin(), vy.end())};
            wall.bbox_hi = {*std::max_element(vx.begin(), vx.end()),
                            *std::max_element(vy.begin(), vy.end())};
        }
        else if (shape == "parser")
        {
            pp_eb.get("implicit_function", wall.implicit_function);
            auto const lo = get_array(pp_eb, "bounding_box_lo");
            auto const hi = get_array(pp_eb, "bounding_box_hi");
            if (lo.size() != 2 || hi.size() != 2 || lo[0] >= hi[0] || lo[1] >= hi[1]) {
                throw std::runtime_error(
                    "eb.bounding_box_lo/hi must be two values each (x y), with lo < hi");
            }
            wall.bbox_lo = {lo[0], lo[1]};
            wall.bbox_hi = {hi[0], hi[1]};
        }
        else
        {
            throw std::runtime_error(
                "eb.shape must be none, elliptical, rectangular, polygon or parser, but is: "
                + shape);
        }

#ifndef AMREX_USE_EB
        throw std::runtime_error("eb.shape requires ImpactX compiled with ImpactX_EB=ON");
#endif
        return wall;
    }

    void
    build_eb (AmrCoreData & amr_data)
    {
        BL_PROFILE("impactx::initialization::build_eb");

        clear_eb(amr_data);

        auto const wall = read_transverse_wall();
        if (!wall) {
            return;
        }

#ifdef AMREX_USE_EB
        amrex::ParmParse const pp_eb("eb");
        amrex::Parser const parser = pp_eb.makeParser(wall->implicit_function, {"x", "y", "z"});
        auto const shop = amrex::EB2::makeShop(amrex::EB2::ParserIF(parser.compile<3>()), parser);

        auto & tp = amr_data.track_particles;
        for (int lev = 0; lev <= amr_data.finestLevel(); ++lev)
        {
            amrex::Geometry const & geom = amr_data.Geom(lev);

            // no coarse EB levels: the multigrid solver coarsens the level set and edge centroids
            amrex::EB2::Build(shop, geom, 0, 0);
            amrex::EB2::IndexSpace * index_space = &amrex::EB2::IndexSpace::top();
            tp.m_eb_index_space[lev] = index_space;

            tp.m_eb_factory[lev] = amrex::makeEBFabFactory(
                index_space, geom, amr_data.boxArray(lev), amr_data.DistributionMap(lev),
                {2, 2, 2}, amrex::EBSupport::full);
        }
#endif
    }

    amrex::Vector<EBFactory const *>
    eb_factories ([[maybe_unused]] AmrCoreData const & amr_data)
    {
        amrex::Vector<EBFactory const *> factories;
#ifdef AMREX_USE_EB
        auto const & eb_factory = amr_data.track_particles.m_eb_factory;
        if (!eb_factory.empty()) {
            for (int lev = 0; lev <= amr_data.finestLevel(); ++lev) {
                factories.push_back(eb_factory.at(lev).get());
            }
        }
#endif
        return factories;
    }

    void
    clear_eb ([[maybe_unused]] AmrCoreData & amr_data)
    {
#ifdef AMREX_USE_EB
        auto & tp = amr_data.track_particles;

        // factories point into the index spaces: release them first
        tp.m_eb_factory.clear();
        for (auto const & [lev, index_space] : tp.m_eb_index_space) {
            amrex::EB2::IndexSpace::erase(index_space);
        }
        tp.m_eb_index_space.clear();
#endif
    }

} // namespace impactx::initialization
