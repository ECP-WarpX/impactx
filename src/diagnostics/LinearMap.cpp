/* Copyright 2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "LinearMap.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>

#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>


namespace impactx::diagnostics
{
namespace
{
    /** Evaluate @c element.transport_map(ref), or apply the
     *  missing-linear-map policy for elements whose @c transport_map is
     *  known not to be implemented.
     *
     *  The branch is chosen at compile time via
     *  @c T_Element::has_linear_transport (provided by the
     *  @c mixin::LinearTransport base):
     *
     *    - Elements with an implemented linear map (the vast majority):
     *      the call is forwarded directly without any @c try/@c catch.
     *      This is the hot path and compiles to a single virtual-free
     *      call with no exception machinery.
     *    - Elements without an implemented linear map: the @p on_missing
     *      policy is applied -- either re-raise as a
     *      @c std::runtime_error naming the element, substitute the
     *      identity with a one-per-type warning, or substitute the
     *      identity silently.
     */
    template <typename T_Element>
    Map6x6
    safe_transport_map (
        T_Element const & element,
        RefPart const & ref,
        std::set<std::string> & warned_types,
        OnMissingLinearMap on_missing
    )
    {
        if constexpr (T_Element::has_linear_transport)
        {
            // Hot path: element has an implemented linear map. Avoid the
            // parameter/set/string work in the "missing" branch entirely
            // so the compiler can inline this to a direct call.
            amrex::ignore_unused(warned_types, on_missing);
            // the transfer map is physics; the name used in the cold path below is not
            return elements::physics_of(element).transport_map(ref);
        }
        else
        {
            // Cold path: only instantiated for element types whose
            // transport_map is known to throw. No try/catch needed.
            std::string const type_name{T_Element::type};
            switch (on_missing)
            {
                case OnMissingLinearMap::Throw:
                {
                    std::string msg = "Undefined linear transport map in lattice for element ";
                    if (element.has_name()) { msg += element.name() + " "; }
                    msg += "of type " + type_name;
                    throw std::runtime_error(msg);
                }
                case OnMissingLinearMap::IdentityWithWarning:
                {
                    if (warned_types.insert(type_name).second)
                    {
                        ablastr::warn_manager::WMRecordWarning(
                            "ImpactX::diagnostics::map_trace",
                            "Element type '" + type_name + "' does not provide a "
                            "linear transport map. It is being treated as the "
                            "identity for this diagnostic. Any linear-optics "
                            "output (Twiss, tune, chromaticity, dispersion) will "
                            "therefore not include this element's contribution.",
                            ablastr::warn_manager::WarnPriority::medium
                        );
                    }
                    return Map6x6::Identity();
                }
                case OnMissingLinearMap::IdentitySilent:
                {
                    // User explicitly opted in to identity fallback;
                    // suppress the warning as they already know.
                    return Map6x6::Identity();
                }
            }
            // All enum values handled above; unreachable.
            return Map6x6::Identity();
        }
    }

    /** Functor type used as the default @c walk_lattice exit hook.
     *
     *  Its call operator does nothing. When used as the default template
     *  argument of @ref walk_lattice, it allows callers that only care
     *  about the cumulative linear map to omit the hook altogether.
     */
    struct NoopElementExit
    {
        template <typename T_Element>
        void operator() (
            T_Element const & /*element*/,
            RefPart const & /*ref*/,
            Map6x6 const & /*M_cum*/
        ) const noexcept
        {}
    };

    /** Traverse the lattice, advancing the reference particle and
     *  composing the cumulative linear transport map, with a
     *  per-element exit hook.
     *
     *  The reference particle passed in is copied internally; the
     *  caller's reference particle is never modified. The cumulative
     *  map is returned, starting from the identity; a fresh warning-
     *  deduplication set is maintained per call.
     *
     *  Uses the same per-element advance / linear-map evaluation
     *  convention as the lattice traversal in
     *  @c src/tracking/common.H , so the resulting linear optics
     *  is consistent with the tracker for edge-sensitive elements
     *  such as @c RFCavity, @c SoftQuad, and @c SoftSol.
     *
     * @tparam F_OnElementExit callable
     *                         @c void(E const&, RefPart const&, Map6x6 const&)
     *                         invoked once per element after all
     *                         slices have been processed; defaults to
     *                         @c NoopElementExit (no-op)
     * @param[in] lattice         lattice to traverse
     * @param[in] ref_part_init   reference particle at the starting s
     *                            (copied internally)
     * @param[in] on_missing      policy for elements without a
     *                            linear transport map
     * @param[in] on_element_exit per-element hook (optional)
     * @return end-to-end 6x6 linear transport map
     */
    template <typename F_OnElementExit = NoopElementExit>
    Map6x6
    walk_lattice (
        Lattice const & lattice,
        RefPart const & ref_part_init,
        OnMissingLinearMap on_missing,
        F_OnElementExit && on_element_exit = NoopElementExit{}
    )
    {
        // private copy of the reference particle
        RefPart ref = ref_part_init;
        Map6x6 M_cum = Map6x6::Identity();

        // warn only once per element type
        std::set<std::string> warned_types;

        // The element loop mirrors the lattice traversal in src/tracking/common.H,
        // without its period loop and collective-effect slices
        for (auto const & element_variant : lattice)
        {
            ref.set_edge();

            elements::visit([&](auto && element)
            {
                using E = std::decay_t<decltype(element)>;

                // unwrap the host-only metadata: reference tracking is physics, not
                // book-keeping. safe_transport_map() keeps the full element, because its
                // diagnostics name the element.
                auto && physics = elements::physics_of(element);

                int const nslice = element.nslice();
                for (int slice_i = 0; slice_i < nslice; ++slice_i)
                {
                    physics(ref);
                    Map6x6 const R_slice = safe_transport_map<E>(
                        element, ref, warned_types, on_missing
                    );
                    M_cum = R_slice * M_cum;
                }

                on_element_exit(element, ref, M_cum);
            }, element_variant);
        }

        return M_cum;
    }
} // namespace

    std::vector<MapTraceEntry>
    map_trace (
        Lattice const & lattice,
        RefPart const & ref_part_init,
        OnMissingLinearMap on_missing
    )
    {
        BL_PROFILE("impactx::diagnostics::map_trace");

        std::vector<MapTraceEntry> trace;
        trace.reserve(lattice.size() + 1u);

        // Leading entry: identity map at the starting s. This makes the
        // returned vector always start at the beginning of the lattice so
        // downstream code (e.g., Twiss propagation) can treat every entry uniformly.
        trace.push_back(MapTraceEntry{
            ref_part_init.s,
            std::string{},
            std::string{"<start>"},
            Map6x6::Identity()
        });

        walk_lattice(
            lattice, ref_part_init, on_missing,
            [&trace](auto const & element, RefPart const & ref_now, Map6x6 const & M_now)
            {
                using E = std::decay_t<decltype(element)>;
                trace.push_back(MapTraceEntry{
                    ref_now.s,
                    element.has_name() ? element.name() : std::string{},
                    std::string{E::type},
                    M_now
                });
            }
        );

        return trace;
    }

    Map6x6
    linear_map (
        Lattice const & lattice,
        RefPart const & ref_part_init,
        OnMissingLinearMap on_missing
    )
    {
        BL_PROFILE("impactx::diagnostics::linear_map");

        // End-to-end map only: rely on walk_lattice's default no-op
        // element-exit hook (NoopElementExit), so no per-element
        // MapTraceEntry (two std::string copies plus a 6x6 matrix per
        // element) is constructed.
        return walk_lattice(lattice, ref_part_init, on_missing);
    }

} // namespace impactx::diagnostics
