/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "elements/transformation/Insert.H"

#include "elements/helper/Accessors.H"

#include <optional>
#include <stdexcept>
#include <utility>


namespace impactx::elements::transformation
{
    Lattice
    insert_element_every_ds (
        Lattice const & lattice,
        amrex::ParticleReal ds,
        elements::KnownElements element
    )
    {
        // algorithm below is so far only implemented for thin elements to insert
        double const new_element_ds = elements::ds(element);  // in meters
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            new_element_ds == 0,
            "insert_element_ever_s: Only thin elements are supported."
        );

        Lattice new_lattice;

        double s = 0.0;  // in meters   // TODO: if we can avoid a global s, we can avoid wasting significant digits for long lattices
        double s_next_insert = ds;  // in meters

        // the tail of an element that was split, still waiting to be placed
        std::optional<elements::ElementHandle> pending;
        Lattice::size_type next = 0;

        while (pending.has_value() || next < lattice.size())
        {
            // take the leftover of the last split first, else the next element
            elements::ElementHandle cur_element = pending.has_value()
                ? *std::exchange(pending, std::nullopt)
                : lattice[next++];

            // check where the current element ends
            double const cur_s_out = s + elements::ds(cur_element);  // in meters

            // case 1: current element is thick and ends after next insert
            if (s_next_insert < cur_s_out)
            {
                double const s_rel_insert = s_next_insert - s;

                if (elements::is_thin(cur_element))
                {
                    throw std::runtime_error("insert_element_ever_s: Thin element cannot be split.");
                }

                // splitting creates two new physical elements, so neither may alias the
                // element that was split: the caller still owns that one
                elements::ElementHandle head = elements::copy_element(cur_element);
                elements::ElementHandle leftover = elements::copy_element(cur_element);

                elements::ds(head, static_cast<amrex::ParticleReal>(s_rel_insert));
                elements::ds(
                    leftover,
                    elements::ds(leftover) - static_cast<amrex::ParticleReal>(s_rel_insert)
                );
                elements::name(leftover, elements::name(leftover) + "_leftover");

                // insert element in between
                new_lattice.push_back(std::move(head));
                new_lattice.emplace_back(element);

                // the tail is carried into the next iteration
                pending = std::move(leftover);

                s += s_rel_insert;
                s_next_insert += ds;
            }
            // case 2: current element ends exactly with next insert
            else if (s_next_insert == cur_s_out) {
                new_lattice.push_back(std::move(cur_element));
                new_lattice.emplace_back(element);

                s = cur_s_out;
                s_next_insert += ds;
            }
            // case 3: current element ends before next insert
            else {
                // thin element or element too thin to slice in ds
                new_lattice.push_back(std::move(cur_element));

                s = cur_s_out;
            }
        }

        return new_lattice;
    }

} // namespace impactx::elements::transformation
