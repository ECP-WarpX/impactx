#!/usr/bin/env python
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Eric G. Stern
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

# Create and return a list of ImpactX elements from a
# Synergia lattice

# This probably has to be called after all the ImpactX initialization with
# the reference particle.

from enum import Enum

import numpy as np

import impactx

from ..extensions.KnownElementsList import _filter_kwargs
from . import synmadx_pybind as synmadx


class Order(Enum):
    linear = 0
    chr = 1
    exact = 2


# set the following switch overrides the model for dipedges
force_linear_dipedges_for_bends = True


ET = synmadx.element_type

# nslice_by_elem_type = {
#     "drift": 1,
#     "sbend": 6,
#     "quadrupole": 8,
#     "sextupole": 1,
#     "octupole": 1,
#     "multipole": 1,
#     "rfcavity": 1,
# }
nslice_by_elem_type = {
    "drift": 1,
    "sbend": 1,
    "quadrupole": 1,
    "sextupole": 1,
    "octupole": 1,
    "multipole": 1,
    "rfcavity": 1,
}


# linear flag if true uses linearized models
def cnv_drift(elem, order):
    ds = elem.get_length()
    nm = elem.get_name()
    ns = nslice_by_elem_type["drift"]
    if order == Order.exact:
        return impactx.elements.ExactDrift(ds, nslice=ns, name=nm)
    elif order == Order.linear:
        return impactx.elements.Drift(ds, nslice=ns, name=nm)
    elif order == Order.chr:
        return impactx.elements.ChrDrift(ds, nslice=ns, name=nm)
    else:
        raise RuntimeError(f"unknown order: {order}")


def cnv_dipedge(elem, order):
    rc = 1 / elem.get_double_attribute("h", 0.0)
    e1 = elem.get_double_attribute("e1", 0.0)
    fint = elem.get_double_attribute("fint", 0.0)
    hgap = elem.get_double_attribute("hgap", 0.0)
    # sneaky, MAD-X dipedge uses HGAP or half-gap while ImpactX uses
    # g for full gap
    tilt = elem.get_double_attribute("tilt", 0.0)
    rot = 180 * tilt / np.pi
    model = " "
    if order == Order.linear:
        model = "linear"
    elif order == Order.exact or order == Order.chr:
        model = "nonlinear"
    ##
    ## ImpactX has the location= argument which can be either 'entry' or 'exit'
    ## but MAD-X doesn't. How do we resolve this?
    ix_elem = impactx.elements.DipEdge(
        psi=e1,
        rc=rc,
        g=hgap * 2,
        K2=fint,
        model=model,
        name=elem.get_name(),
        rotation=rot,
    )
    return ix_elem


# linear flag if true uses linearized ImpactX models
def cnv_sbend(elem, order):
    bendangle = elem.get_bend_angle()
    length = elem.get_length()
    radius_of_curvature = length / bendangle
    nm = elem.get_name()
    k1 = elem.get_double_attribute("k1", 0.0)
    k1s = elem.get_double_attribute("k1s", 0.0)
    k2 = elem.get_double_attribute("k2", 0.0)
    e1 = elem.get_double_attribute("e1", 0.0)
    e2 = elem.get_double_attribute("e2", 0.0)
    fint = elem.get_double_attribute("fint", 0.0)
    hgap = elem.get_double_attribute("hgap", 0.0)
    tilt = elem.get_double_attribute("tilt", 0.0)
    rot = 180 * tilt / np.pi

    us_dipedge = None
    ds_dipedge = None
    cf = (k1 != 0.0) or (k2 != 0.0) or (k1s != 0.0)

    # model for dipedges
    if force_linear_dipedges_for_bends or order == Order.linear:
        de_model = "linear"
    else:  # Order.exact, Order.chr
        de_model = "nonlinear"

    if e1 != 0.0:
        us_dipedge = impactx.elements.DipEdge(
            e1,
            radius_of_curvature,
            2 * hgap,
            K2=fint,
            location="entry",
            model=de_model,
            name=nm + "_usedge",
            rotation=rot,
        )
    if e2 != 0.0:
        ds_dipedge = impactx.elements.DipEdge(
            e2,
            radius_of_curvature,
            2 * hgap,
            K2=fint,
            location="exit",
            model=de_model,
            name=nm + "_dsedge",
            rotation=rot,
        )

    if elem.get_double_attribute("h1", 0.0) != 0.0:
        print("h1 attribute not supported for sbend in ImpactX")
    if elem.get_double_attribute("h2", 0.0) != 0.0:
        print("h2 attribute not supported for sbend in ImpactX")

    ns = nslice_by_elem_type["sbend"]
    ds = length

    if not cf:
        # normal bend, no higher order moments
        # What kind of screwy program accepts angles in degrees?

        if order == Order.linear:
            main_bend_elem = impactx.elements.Sbend(
                ds, radius_of_curvature, name=nm, nslice=ns, rotation=rot
            )
        else:
            phi = bendangle * 180 / np.pi
            main_bend_elem = impactx.elements.ExactSbend(
                ds, phi, nslice=ns, rotation=rot, name=nm
            )
    else:
        if order == Order.linear:
            # only do first order bend
            main_bend_elem = impactx.elements.CFbend(
                ds, radius_of_curvature, k1, nslice=ns, name=nm, rotation=rot
            )
        else:
            # CF bend
            if k2 == 0.0:
                knormal = [1 / radius_of_curvature, k1]
                kskew = [0.0, k1s]
            else:
                knormal = [1 / radius_of_curvature, k1, k2]
                kskew = [0.0, k1s, 0.0]

            main_bend_elem = impactx.elements.ExactCFbend(
                ds=length,
                k_normal=knormal,
                k_skew=kskew,
                int_order=6,
                mapsteps=6,
                nslice=ns,
                name=nm,
                rotation=rot,
            )

    # collect the pieces
    ixelem = []
    if us_dipedge:
        ixelem.append(us_dipedge)
    ixelem.append(main_bend_elem)
    if ds_dipedge:
        ixelem.append(ds_dipedge)

    return ixelem


def cnv_rbend(elem, order):
    # the RBEND is converted to an SBEND with dipedge elements
    # before an after to get parallel faces

    bendangle = elem.get_bend_angle()
    straight_length = elem.get_length()
    # MAD-X RBEND langths is the straight length.
    radius_of_curvature = straight_length / (2 * np.sin(0.5 * bendangle))
    length = radius_of_curvature * bendangle

    nm = elem.get_name()
    k1 = elem.get_double_attribute("k1", 0.0)
    k1s = elem.get_double_attribute("k1s", 0.0)
    k2 = elem.get_double_attribute("k2", 0.0)
    e1 = elem.get_double_attribute("e1", 0.0)
    e2 = elem.get_double_attribute("e2", 0.0)
    fint = elem.get_double_attribute("fint", 0.0)
    hgap = elem.get_double_attribute("hgap", 0.0)
    tilt = elem.get_double_attribute("tilt", 0.0)
    rot = 180 * tilt / np.pi
    us_dipedge = None
    ds_dipedge = None
    cf = (k1 != 0.0) or (k2 != 0.0) or (k1s != 0.0)

    # for RBENDS, the dipedge angles are relative to bendangle/2 and
    # have opposite sense depending on the sign of the bend.
    if bendangle > 0.0:
        e1 = e1 + bendangle / 2
        e2 = e2 + bendangle / 2
    else:
        e1 = -e1 - bendangle / 2
        e2 = -e2 - bendangle / 2

    # model for dipedges
    if force_linear_dipedges_for_bends or order == Order.linear:
        de_model = "linear"
    else:  # Order.exact, Order.chr
        de_model = "nonlinear"

    if e1 != 0.0:
        us_dipedge = impactx.elements.DipEdge(
            e1,
            radius_of_curvature,
            2 * hgap,
            K2=fint,
            model=de_model,
            name=nm + "_usedge",
            rotation=rot,
        )

    if e2 != 0.0:
        ds_dipedge = impactx.elements.DipEdge(
            e2,
            radius_of_curvature,
            2 * hgap,
            K2=fint,
            model=de_model,
            name=nm + "_dsedge",
            rotation=rot,
        )

    if elem.get_double_attribute("h1", 0.0) != 0.0:
        print("h1 attribute not supported for sbend in ImpactX")
    if elem.get_double_attribute("h2", 0.0) != 0.0:
        print("h2 attribute not supported for sbend in ImpactX")

    ns = nslice_by_elem_type["sbend"]
    ds = length

    if not cf:
        # normal bend, no higher order moments
        # What kind of screwy program accepts angles in degrees?

        if order == Order.linear:
            main_bend_elem = impactx.elements.Sbend(
                ds, radius_of_curvature, name=nm, nslice=ns, rotation=rot
            )
        else:
            phi = bendangle * 180 / np.pi
            main_bend_elem = impactx.elements.ExactSbend(
                ds, phi, nslice=ns, name=nm, rotation=rot
            )

    else:
        # CF bend

        if order == Order.linear:
            # only do linear bend
            main_bend_elem = impactx.elements.CFbend(
                ds, radius_of_curvature, k1, nslice=ns, name=nm
            )
        else:
            # full CF bend
            if k2 == 0.0:
                knormal = [1 / radius_of_curvature, k1]
                kskew = [0.0, k1s]
            else:
                knormal = [1 / radius_of_curvature, k1, k2]
                kskew = [0.0, k1s, 0.0]

            main_bend_elem = impactx.elements.ExactCFbend(
                ds=length,
                k_normal=knormal,
                k_skew=kskew,
                int_order=2,
                nslice=ns,
                name=nm,
                rotation=rot,
            )

    # collect the pieces
    ixelem = []
    if us_dipedge:
        ixelem.append(us_dipedge)
    ixelem.append(main_bend_elem)
    if ds_dipedge:
        ixelem.append(ds_dipedge)

    return ixelem


# linear flag if true uses linearized models
def cnv_quadrupole(elem, order):
    ds = elem.get_length()
    k1 = elem.get_double_attribute("k1", 0.0)
    nm = elem.get_name()
    tilt = elem.get_double_attribute("tilt", 0.0)
    rot = 180 * tilt / np.pi
    if order == Order.linear:
        ix_elem = impactx.elements.Quad(
            ds, k1, nslice=nslice_by_elem_type["quadrupole"], name=nm, rotation=rot
        )
    elif order == Order.chr:
        ix_elem = impactx.elements.ChrQuad(
            ds, k1, nslice=nslice_by_elem_type["quadrupole"], name=nm, rotation=rot
        )
    elif order == Order.exact:
        ix_elem = impactx.elements.ExactQuad(
            ds,
            k1,
            nslice=nslice_by_elem_type["quadrupole"],
            mapsteps=6,
            int_order=6,
            name=nm,
            rotation=rot,
        )
    else:
        raise RuntimeError(f"error, unknown order: {order}")

    return ix_elem


def cnv_multipole(elem, order):
    # ImpactX multipole elements only have a single order each so
    # we have to peel each order and possibly create multiple elements to get the
    # full description.
    nm = elem.get_name()
    knlvect = elem.get_vector_attribute("knl", [])
    kslvect = elem.get_vector_attribute("ksl", [])
    nlen = len(knlvect)
    slen = len(kslvect)
    maxlen = max(nlen, slen)
    tilt = elem.get_double_attribute("tilt", 0.0)
    rot = 180 * tilt / np.pi
    if order == Order.linear:
        # truncate at order 1
        maxlen = max(maxlen, 2)  # 2 means dipole+quadrupole
    norm_mom = np.zeros(maxlen)
    skew_mom = np.zeros(maxlen)
    norm_mom[:nlen] = knlvect
    skew_mom[:slen] = kslvect
    elem_list = []
    for i in range(maxlen):
        # is there anything this order?
        if norm_mom[i] == 0.0 and skew_mom[i] == 0.0:
            pass
        elem_list.append(
            impactx.elements.Multipole(
                i + 1, norm_mom[i], skew_mom[i], name=f"{nm}_{order}", rotation=rot
            )
        )
    return elem_list


def cnv_rfcavity(elem, refpart, order):
    # I don't have a linear model of the RF cavity itself, but if there are
    # drifts surrounding it, they should be linear if the linear model is
    # requested
    mp = refpart.get_mass()
    # The cavity voltage is in relationship to the reference particle mass
    # The MAD-X convention is for the voltage to be in MV
    rfvolt = elem.get_double_attribute("volt", 0.0) / 1000.0  # convert to GV
    freq = elem.get_double_attribute("freq", 0.0) * 1.0e6  # get the freq in Hz
    phase = elem.get_double_attribute("lag", 0.0) * 360.0 - 90.0
    # if cavity length > 0, create two drifts to sandwich it
    L = elem.get_length()

    # if cavity length > 0, create two drifts to sandwich it
    RFelem = impactx.elements.ShortRF(rfvolt / mp, freq, phase, name=elem.get_name())
    if L == 0.0:
        # exactly 0 should be allowable for == comparison
        return RFelem
    else:
        halfL = L / 2
        if order == Order.linear:
            d1 = impactx.elements.Drift(halfL, nslice=1, name=f"{elem.get_name()}_U")
            d2 = impactx.elements.Drift(halfL, nslice=1, name=f"{elem.get_name()}_D")
        else:
            d1 = impactx.elements.ExactDrift(
                halfL, nslice=1, name=f"{elem.get_name()}_U"
            )
            d2 = impactx.elements.ExactDrift(
                halfL, nslice=1, name=f"{elem.get_name()}_D"
            )

        RFunit = [d1, RFelem, d2]
    return RFunit


def cnv_kicker(elem, order):
    L = elem.get_length()
    ename = elem.get_name()

    xkick = 0.0
    ykick = 0.0
    if elem.get_type() == ET.hkicker:
        if elem.has_double_attribute("kick"):
            xkick = elem.get_double_attribute("kick")
        elif elem.has_double_attribute("hkick"):
            # some lattice files use hkick attribute instead of kick
            xkick = elem.get_double_attribute("hkick")

    elif elem.get_type() == ET.vkicker:
        if elem.has_double_attribute("kick"):
            ykick = elem.get_double_attribute("kick")
        elif elem.has_double_attribute("vkick"):
            # some lattice files use vkick attribute instead of kick
            ykick = elem.get_double_attribute("vkick")

    elif elem.get_type() == ET.kicker:
        if elem.has_double_attribute("hkick"):
            xkick = elem.get_double_attribute("hkick")
        if elem.has_double_attribute("vkick"):
            ykick = elem.get_double_attribute("vkick")

    kickelem = impactx.elements.Kicker(xkick=xkick, ykick=ykick, name=ename)

    # if length is exactly 0, we're done.
    if L == 0.0:
        return kickelem
    # if length is >0, create drifts on both sides of the
    # thin kick

    if order == Order.linear:
        d1 = impactx.elements.Drift(L / 2, nslice=1, name=f"{ename}_U")
        d2 = impactx.elements.Drift(L / 2, nslice=1, name=f"{ename}_D")
    else:
        d1 = impactx.elements.ExactDrift(L / 2, nslice=1, name=f"{ename}_U")
        d2 = impactx.elements.ExactDrift(L / 2, nslice=1, name=f"{ename}_D")

    kickunit = [d1, kickelem, d2]
    return kickunit


# Hi Rob,

# Here's the full expression for  converting MAD8 multipoles to MADX.

# For MAD8 multipole KnL=str and Tn=angle where n could be 1,2,3,4, ....

# The MADX normal component is COS((n+1)*angle)
# The MADX skew component is -SIN((n+1)*angle)

# Yes, that is a negative sign only on the skew component because the
# derivation contains a rotation around the negative angle.

# MAD8 multipole
# MULTIPOLE, K2L=str, T2=angle

# Converts to MADX statement
# MULTIPOLE, KNL={0, 0, str*cos(3*angle)}, KSL={0, 0, -str*sin(3*angle))};

# Your usual multipole elements that specify an angle in the Recycler
# latticelike the one below  have written out values for angles, so for
# instance T2=pi/6 which is the normal skew angle for sextupoles. The normal
# component vanishes and the skew component receives negative the strength.

# MP114BS: MULTIPOLE, TYPE=RGF021_S_BODY, K1L=0.16E-4, T1=0.785398163397,
# K2L=&
# 0.808E-3*DIR, T2=0.523598775598, K3L=0.064894, T3=0.392699081699, K4L=&
# 15.028637*DIR, T4=0.314159265359

# Converts to
# MP114BS: MULTIPOLE, TYPE=RGF021_S_BODY, KSL={-0.16E-4,  -0.808E-3*DIR,
# -0.064894, -15.028637*DIR};


#     Eric


def cnv_sextupole(elem, order):
    L = elem.get_length()
    k2 = elem.get_double_attribute("k2", 0.0)
    # The Booster lattice includes elements with the tilt attribute
    tilt = elem.get_double_attribute("tilt", 0.0)
    # The tilt attribute is encoded into the multipole moments
    k2n = k2 * np.cos(3 * tilt)
    k2s = -k2 * np.sin(3 * tilt)
    nm = elem.get_name()
    knorm = np.array([0, 0, k2n])
    kskew = np.array([0, 0, k2s])

    if order == Order.linear:
        # There is no linear Sextupole. Use a linear drift
        sxelem = impactx.elements.Drift(L, name=nm)
    else:
        sxelem = impactx.elements.ExactMultipole(
            ds=L,
            k_normal=knorm,
            k_skew=kskew,
            int_order=4,
            nslice=nslice_by_elem_type["sextupole"],
            name=nm,
        )
    return sxelem


# octupole follows similar logic as sextupole
def cnv_octupole(elem, order):
    L = elem.get_length()
    k3 = elem.get_double_attribute("k3", 0.0)
    # The Booster lattice includes elements with the tilt attribute
    # The tilt value is encoded into the multipole moment
    tilt = elem.get_double_attribute("tilt", 0.0)
    k3n = k3 * np.cos(4 * tilt)
    k3s = -k3 * np.sin(4 * tilt)
    nm = elem.get_name()
    knorm = np.array([0, 0, 0, k3n])
    kskew = np.array([0, 0, 0, k3s])

    if order == Order.linear:
        # There is no linear Octupole. Use a linear drift
        ocelem = impactx.elements.Drift(L, name=nm)
    else:
        ocelem = impactx.elements.ExactMultipole(
            ds=L,
            k_normal=knorm,
            k_skew=kskew,
            int_order=4,
            nslice=nslice_by_elem_type["octupole"],
            name=nm,
        )

    return ocelem


def cnv_monitor(elem):
    fname = elem.get_string_attribute("file", "monitor")
    ix_elem = impactx.elements.BeamMonitor(name=fname, backend="h5")
    return ix_elem


def cnv_marker(elem):
    ix_elem = impactx.elements.Drift(ds=0.0, name=elem.get_name())
    return ix_elem


def syn2_to_impactx(lattice, init_monitor=True, final_monitor=True, order=Order.exact):
    # lattice must have a reference particle
    try:
        refpart = lattice.get_reference_particle()
    except RuntimeError:
        print("cannot get reference particle.")
        return None
    # print("using order:", order)

    impactx_lattice = []

    # We may define the monitor element if needed
    monitor = None
    if init_monitor:
        if not monitor:
            monitor = impactx.elements.BeamMonitor("monitor", backend="h5")
        impactx_lattice.append(monitor)

    # peel elements from the synergia lattice, converting to ImpactX elements
    for elem in lattice.get_elements():
        etype = elem.get_type()

        if etype == ET.drift:
            impactx_lattice.append(cnv_drift(elem, order))
        elif etype == ET.sbend:
            bndelem = cnv_sbend(elem, order)
            if isinstance(bndelem, list):
                impactx_lattice.extend(bndelem)
            else:
                impactx_lattice.append(bndelem)
        elif etype == ET.rbend:
            bndelem = cnv_rbend(elem, order)
            if isinstance(bndelem, list):
                impactx_lattice.extend(bndelem)
            else:
                impactx_lattice.append(bndelem)
        elif etype == ET.quadrupole:
            impactx_lattice.append(cnv_quadrupole(elem, order))
        elif etype == ET.sextupole:
            sxelem = cnv_sextupole(elem, order)
            if isinstance(sxelem, list):
                impactx_lattice.extend(sxelem)
            else:
                impactx_lattice.append(sxelem)
        elif etype == ET.octupole:
            ocelem = cnv_octupole(elem, order)
            if isinstance(ocelem, list):
                impactx_lattice.extend(ocelem)
            else:
                impactx_lattice.append(ocelem)
        elif etype == ET.hkicker or etype == ET.vkicker or etype == ET.kicker:
            # both H and V kickers handled by same routine since
            # ImpactX has only one type of kicker element
            impactx_lattice.append(cnv_kicker(elem, order))
        elif etype == ET.nllens:
            raise RuntimeError("nllens not yet implemented")
            # impactx_lattice.append(cnv_nllens(elem))
        elif etype == ET.dipedge:
            impactx_lattice.append(cnv_dipedge(elem, order))
        elif etype == ET.rfcavity:
            rfelem = cnv_rfcavity(elem, refpart, order)
            if isinstance(rfelem, list):
                impactx_lattice.extend(rfelem)
            else:
                impactx_lattice.append(rfelem)
        elif etype == ET.multipole:
            mpelem = cnv_multipole(elem, order)
            if isinstance(mpelem, list):
                impactx_lattice.extend(mpelem)
            else:
                impactx_lattice.append(mpelem)
        elif etype == ET.monitor:
            ix_elem = cnv_monitor(elem)
            impactx_lattice.append(ix_elem)
        elif etype == ET.marker:
            ix_elem = cnv_marker(elem)
            impactx_lattice.append(ix_elem)
        else:
            print("warning: unsupported element: ", etype)

    if final_monitor:
        # define monitor if it hasn't already been defined
        if not monitor:
            monitor = impactx.elements.BeamMonitor("monitor", backend="h5")
        impactx_lattice.append(monitor)

    return impactx_lattice


# given an ImpactX lattice as a sequence of elements, unroll the lattice
# and return a python program that reccreates the lattice


def unroll_impactx_lattice(lattice):
    # define the output lattice as a list

    output_lattice = "[\n"
    for elem in lattice:
        output_elem = ""
        edict = elem.to_dict()
        etype = edict["type"]
        output_elem = f"impactx.elements.{etype}("
        firstparm = True

        # for ExactSbend the member phi which is the angle has been converted into radians but
        # I need to convert it back to degrees which is what the contructor needs.
        if etype == "ExactSbend":
            edict["phi"] = edict["phi"] * 180 / np.pi

        # drop what the constructor does not take (``type``; ``ds`` of thin elements)
        edict = _filter_kwargs(edict)

        # skipping BeamMonitors for now. They seem to cause trouble
        if etype == "BeamMonitor":
            continue

        for pname in edict:
            # the name parameter is a string and must be enclosed in quotes,
            # unless the element is unnamed (None).
            # Also location and model parameters on the DipEdge element
            if not firstparm:
                output_elem = output_elem + ", "
            if pname == "name" and edict.get(pname) is None:
                output_elem = output_elem + f"{pname}=None"
            elif pname == "name" or pname == "location" or pname == "model":
                output_elem = output_elem + f'{pname}="{edict.get(pname)}"'
            else:
                output_elem = output_elem + f"{pname}={edict.get(pname)}"
            firstparm = False
        # close out this element
        output_elem = output_elem + ")"

        output_lattice = output_lattice + output_elem + ",\n"

    # close out the lattice
    output_lattice = output_lattice + "]"

    return output_lattice


__misc_txt = """
    beam, particle=proton, energy=0.8+pmass;
    d: drift, l=1.0;
    b: sbend, angle=pi/12.0;
    rfc: rfcavity, l=2.4, volt=0.02, freq=53.0;
    q1: quadrupole, l=0.5, k1=1.0/(0.5*7.0);
    sx1: sextupole, l=0.2, k2=0.05;
    bb: rbend, l=1.5, angle=(2*pi)/96, k1=0.004, k2=0.0005;
    m1: marker;
    m2: marker;
    mon1: monitor;
    sbup: sbend, l=2.0, angle=pi/12, tilt=-pi/2;
    fmag: rbend, l=2.4, angle=pi/15, k1=0.02, k2=-0.3, tilt=pi/6;
    misc: line=(d, m1, mon1, b, rfc, q1, sx1, m2, bb, sbup, fmag);
"""

__simple_booster_txt = """
! Simplified Booster lattice

// From JFO 2022-12-08
// The simplified booster lattice is trivial. I have no madx lattice (you could make one
// very easily) - I just use pyorbit classes directly to instantiate a basic cell;  it is then
// replicated 24 times. I took the bending magnets lengths and
// strengths directly from the official MADX  lattice file.

// The basic cell is

// d1 fmag d2 dmag d3

// d1, d2, d3 : drifts of lengths 0.6 0.5 and 3.0 m
// fmag:  focusing      bend   L = 2.889612 m
// dmag   defocusing bend   L = 2.889612 m

// total cell length: 19.758 m
// total ring  length = 24*19.758 = 474.20 m

// The length, focusing strengths and curvature radius of the
// magnets are as in the booster MADX file.

// If you entered 1 cell correctly, you should get the periodic solution:
// bx = 33.86 m ax = 0
// by = 5.39m    ay =0
// For 24 cells, the raw tunes are nux = 7.017 and nuy = 6.675.  You will need to tweak the nominal focusing strengths a bit to avoid resonances.


//--------- Nominal Gradient Magnet Definitions

// EGS
// The apparent cell structure is actually:

// D1, l=0.6;
// FMAGU01;
// D2, l=0.5;
// DMAGU01;
// D3, l=6.0;
// DMAGD01;
// D4, l=0.5;
// FMAGD01;
// DR, l=0.6


! Expand structure to include corrector packages:

! Normal short straight

! drift 0.6 (end of previous cell)
! drift 0.176
! corrector package (short) l=0.168
! drift 0.256
! end of short straight

! non-RF cavity long straight
!
! drift 5.581
! correction package (long) l=0.168
! drift 0.251
!

! long straight with RF cavity
!
! drift 0.21
! rfcavity drift-cavity-drift 2.35
! drift 0.12
! rfcavity drift-cavity-drift 2.35
! 4 drifts total 0.551
! corrector packages (long) 0.168
! drift 0.251

! Corrector package:
!
! HLxx, HKICKER, l=0.024
! VLxx, VKICKER, l=0.024
! Q{S|L}xx QUADRUPOLE, l=0.024 (normal)
! MULTIPOLE Q{S|L}ERR, l=0
! MONITOR, l=0.024
! QS{S|L}xx, QUADRUPOLE, l=0.024 (skew)
! MULTIPOLE QS{S|L}ERR, l=0
! SEXTUPOLE SX{L|S}, l=0.024 (normal)
! SEXTUPOLE SS{L|S}, l=0.024 (skew)
!

ke1 = 0.8;  !800 MeV kinetic energy at injection

rhof  :=  40.847086;   !  bending radius of focusing magnet
rhod  :=  48.034101;   !  bending radius of defocusing magnet

blength :=     2.889612;    !  arc length for both F and D magnets
blengthf :=    2.889009499; !  physical length (straight length) for F magnet
blengthd :=    2.889176299; !  physical length (straight length) for D magnet


!
! The quad field for the gradient magnet is contained in file " qsdqsf.dat" to be read in before this file !
!
! read from file at time step = 7
qsd := -57.38855012e-3;
qsf := 54.10921561e-3;


! These ssd and ssf strengths come from fitting to 01 Dec 2015 chromaticity data
! and predicts chromaticity much better than using the Drozhdin et al measurements above

ssd :=  -0.04381647074 + ke1*(0.009150934932+ ke1*(-0.0023900895  + ke1*(0.000318068028 -  ke1* 1.6353205e-05)));

ssf :=  -0.006384940088 + ke1*(0.01967542848 + ke1*( -0.006776746 + ke1*(0.00091367565 - ke1* 4.293705e-05)));

 !
 ! Gradient magnets defined by their physical length aith their bend angle
 ! being defined by the arc length/radius of curvature

!FMAG: RBEND,  L = blengthf  , ANGLE = blength/rhof, K1 = qsf  , K2 = ssf;
FMAG: SBEND,  L = blength  , ANGLE = blength/rhof, e1=blength/(2*rhof), e2=blength/(2*rhof), K1 = qsf  , K2 = ssf, type=fmag;
DMAG: SBEND,  L = blength  , ANGLE = blength/rhod, e1=blength/(2*rhod), e2=blength/(2*rhod), K1 = qsd  , K2 = ssd, type=dmag;

! drifts in the short straight section
mins: drift, l=0.5, type=shortstraight;
sc: drift, l=0.6, type=shortstraight;
sa: drift, l=0.176, type=shortstraight;
sb: drift, l=0.256, type=shortstraight;

! drifts in the long straight section
dlong: drift, l=5.581, type=longstraight; ! for no-RF cavity cells
drifta: drift, l=0.21, type=longstraight; ! start of RF cavity cell
driftb: drift, l=0.12, type=longstraight;
dmidls: drift, l=0.551, type=longstraight; ! (drift in the middle of the longstraight)
drifte: drift, l=0.251, type=longstraight; ! end of RF cavity cell

drrf: drift, l=2.35/2, type=rfaperture;
rfc: rfcavity, l=0, harmon=84, volt=0, lag=0, type=rfaperture;

! short corrector package
hsxx: drift, l=0.024, type=shortstraight;
vsxx: drift, l=0.024, type=shortstraight;
qsxx: quadrupole, k1=0.0, l=0.024, type=shortstraight;
! multipole not included
bpms: drift, l=0.024, type=shortstraight;
qssxx: drift, l=0.024, type=shortstraight; // skew quad
! multipole not included
sxsxx: sextupole, k2=0.0, l=0.024, type=shortstraight;
sssxx: drift, l=0.024, type=shortstraight; // skew sextupole

cpshort: line=(hsxx, vsxx, qsxx, bpms, qssxx, sxsxx, sssxx);

! long corrector package
hlxx: drift, l=0.024, type=longstraight;
vlxx: drift, l=0.024, type=longstraight;;
qlxx: quadrupole, k1=0.0, l=0.024, type=shortstraight;
! multipole not included
bpml: drift, l=0.024, type=shortstraight;
qlsxx: drift, l=0.024, type=longstraight; // skew quad
! multipole not included
sxlxx: sextupole, k2=0.0, l=0.024, type=longstraight;
sslxx: drift, l=0.024, type=longstraight; // skew sextupole

cplong: line=(hlxx, vlxx, qlxx, bpms, qlsxx, sxlxx, sslxx);

!!!!!!!!!!!!!!!!   beginning of ring definition

fmagu01: fmag;
fmagd01: fmag;
dmagu01: dmag;
dmagd01: dmag;

cell01 : line = (sa, cpshort, sb, fmagu01, mins, dmagu01, dlong, cplong, drifte, dmagd01, mins, fmagd01, sc);

fmagu02: fmag;
fmagd02: fmag;
dmagu02: dmag;
dmagd02: dmag;

cell02 : line = (sa, cpshort, sb, fmagu02, mins, dmagu02, dlong, cplong, drifte, dmagd02, mins, fmagd02, sc);

fmagu03: fmag;
fmagd03: fmag;
dmagu03: dmag;
dmagd03: dmag;

cell03 : line = (sa, cpshort, sb, fmagu03, mins, dmagu03, dlong, cplong, drifte, dmagd03, mins, fmagd03, sc);

fmagu04: fmag;
fmagd04: fmag;
dmagu04: dmag;
dmagd04: dmag;

cell04 : line = (sa, cpshort, sb, fmagu04, mins, dmagu04, dlong, cplong, drifte, dmagd04, mins, fmagd04, sc);

fmagu05: fmag;
fmagd05: fmag;
dmagu05: dmag;
dmagd05: dmag;

cell05 : line = (sa, cpshort, sb, fmagu05, mins, dmagu05, dlong, cplong, drifte, dmagd05, mins, fmagd05, sc);

fmagu06: fmag;
fmagd06: fmag;
dmagu06: dmag;
dmagd06: dmag;

cell06 : line = (sa, cpshort, sb, fmagu06, mins, dmagu06, dlong, cplong, drifte, dmagd06, mins, fmagd06, sc);

fmagu07: fmag;
fmagd07: fmag;
dmagu07: dmag;
dmagd07: dmag;

cell07 : line = (sa, cpshort, sb, fmagu07, mins, dmagu07, dlong, cplong, drifte, dmagd07, mins, fmagd07, sc);

fmagu08: fmag;
fmagd08: fmag;
dmagu08: dmag;
dmagd08: dmag;

cell08 : line = (sa, cpshort, sb, fmagu08, mins, dmagu08, dlong, cplong, drifte, dmagd08, mins, fmagd08, sc);

fmagu09: fmag;
fmagd09: fmag;
dmagu09: dmag;
dmagd09: dmag;

cell09 : line = (sa, cpshort, sb, fmagu09, mins, dmagu09, dlong, cplong, drifte, dmagd09, mins, fmagd09, sc);

fmagu10: fmag;
fmagd10: fmag;
dmagu10: dmag;
dmagd10: dmag;

cell10 : line = (sa, cpshort, sb, fmagu10, mins, dmagu10, dlong, cplong, drifte, dmagd10, mins, fmagd10, sc);

fmagu11: fmag;
fmagd11: fmag;
dmagu11: dmag;
dmagd11: dmag;

cell11 : line = (sa, cpshort, sb, fmagu11, mins, dmagu11, dlong, cplong, drifte, dmagd11, mins, fmagd11, sc);

fmagu12: fmag;
fmagd12: fmag;
dmagu12: dmag;
dmagd12: dmag;

cell12 : line = (sa, cpshort, sb, fmagu12, mins, dmagu12, dlong, cplong, drifte, dmagd12, mins, fmagd12, sc);

fmagu13: fmag;
fmagd13: fmag;
dmagu13: dmag;
dmagd13: dmag;

cell13 : line = (sa, cpshort, sb, fmagu13, mins, dmagu13, dlong, cplong, drifte, dmagd13, mins, fmagd13, sc);

fmagu14: fmag;
fmagd14: fmag;
dmagu14: dmag;
dmagd14: dmag;
rf01: line=(drrf, rfc, drrf);
rf02: line=(drrf, rfc, drrf);

cell14: line = (sa, cpshort, sb, fmagu14, mins, dmagu14, drifta, rf01, driftb, rf02, dmidls, cplong, drifte, dmagd14, mins, fmagd14, sc);


fmagu15: fmag;
fmagd15: fmag;
dmagu15: dmag;
dmagd15: dmag;
rf03: line=(drrf, rfc, drrf);
rf04: line=(drrf, rfc, drrf);

cell15: line = (sa, cpshort, sb, fmagu15, mins, dmagu15, drifta, rf03, driftb, rf04, dmidls, cplong, drifte, dmagd15, mins, fmagd15, sc);

fmagu16: fmag;
fmagd16: fmag;
dmagu16: dmag;
dmagd16: dmag;
rf05: line=(drrf, rfc, drrf);
rf06: line=(drrf, rfc, drrf);

cell16: line = (sa, cpshort, sb, fmagu16, mins, dmagu16, drifta, rf05, driftb, rf06, dmidls, cplong, drifte, dmagd16, mins, fmagd16, sc);

fmagu17: fmag;
fmagd17: fmag;
dmagu17: dmag;
dmagd17: dmag;
rf07: line=(drrf, rfc, drrf);
rf08: line=(drrf, rfc, drrf);

cell17: line = (sa, cpshort, sb, fmagu17, mins, dmagu17, drifta, rf07, driftb, rf08, dmidls, cplong, drifte, dmagd17, mins, fmagd17, sc);

fmagu18: fmag;
fmagd18: fmag;
dmagu18: dmag;
dmagd18: dmag;
rf09: line=(drrf, rfc, drrf);
rf10: line=(drrf, rfc, drrf);

cell18: line = (sa, cpshort, sb, fmagu18, mins, dmagu18, drifta, rf09, driftb, rf10, dmidls, cplong, drifte, dmagd18, mins, fmagd18, sc);

fmagu19: fmag;
fmagd19: fmag;
dmagu19: dmag;
dmagd19: dmag;
rf11: line=(drrf, rfc, drrf);
rf12: line=(drrf, rfc, drrf);

cell19: line = (sa, cpshort, sb, fmagu19, mins, dmagu19, drifta, rf11, driftb, rf12, dmidls, cplong, drifte, dmagd19, mins, fmagd19, sc);

fmagu20: fmag;
fmagd20: fmag;
dmagu20: dmag;
dmagd20: dmag;
rf13: line=(drrf, rfc, drrf);
rf14: line=(drrf, rfc, drrf);

cell20: line = (sa, cpshort, sb, fmagu20, mins, dmagu20, drifta, rf13, driftb, rf14, dmidls, cplong, drifte, dmagd20, mins, fmagd20, sc);

fmagu21: fmag;
fmagd21: fmag;
dmagu21: dmag;
dmagd21: dmag;
rf15: line=(drrf, rfc, drrf);
rf16: line=(drrf, rfc, drrf);

cell21: line = (sa, cpshort, sb, fmagu21, mins, dmagu21, drifta, rf15, driftb, rf16, dmidls, cplong, drifte, dmagd21, mins, fmagd21, sc);

fmagu22: fmag;
fmagd22: fmag;
dmagu22: dmag;
dmagd22: dmag;
rf17: line=(drrf, rfc, drrf);
rf18: line=(drrf, rfc, drrf);

cell22: line = (sa, cpshort, sb, fmagu22, mins, dmagu22, drifta, rf17, driftb, rf18, dmidls, cplong, drifte, dmagd22, mins, fmagd22, sc);

fmagu23: fmag;
fmagd23: fmag;
dmagu23: dmag;
dmagd23: dmag;
rf19: line=(drrf, rfc, drrf);
rf20: line=(drrf, rfc, drrf);

cell23: line = (sa, cpshort, sb, fmagu23, mins, dmagu23, drifta, rf19, driftb, rf20, dmidls, cplong, drifte, dmagd23, mins, fmagd23, sc);

fmagu24: fmag;
fmagd24: fmag;
dmagu24: dmag;
dmagd24: dmag;
rf21: line=(drrf, rfc, drrf);
rf22: line=(drrf, rfc, drrf);

cell24: line = (sa, cpshort, sb, fmagu24, mins, dmagu24, drifta, rf21, driftb, rf22, dmidls, cplong, drifte, dmagd24, mins, fmagd24, sc);

booster: line=(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24);

beam, particle=proton, energy=0.8+pmass;
!beam, particle=proton, energy=1.7386206689709858;
!beam, particle=proton, energy=0.8003486229709857+pmass;
"""

__iota_txt = """
kqa1r := kq01;
kq01 = -7.72652301;
kqa2r := kq02;
kq02 = 12.28401222;
kqa3r := kq03;
kq03 = -12.43016989;
kqa4r := kq04;
kq04 = 20.16074347;
kqb1r := kq05;
kq05 = -10.24365752;
kqb2r := kq06;
kq06 = 15.12808788;
kqb3r := kq07;
kq07 = -6.92311681;
kqb4r := kq08;
kq08 = -6.90057605;
kqb5r := kq09;
kq09 = 13.50655178;
kqb6r := kq10;
kq10 = -11.91343344;
kqc1r := kq11;
kq11 = -13.51948869;
kqc2r := kq12;
kq12 = 12.0339278;
kqc3r := kq13;
kq13 = -13.56878135;
kqd1r := kq14;
kq14 = -7.97007816;
kqd2r := kq15;
kq15 = 5.92322639;
kqd3r := kq16;
kq16 = -6.32915747;
kqd4r := kq17;
kq17 = 5.1636516;
kqe1r := kq18;
kq18 = -4.69712477;
kqe2r := kq19;
kq19 = 7.0326898;
kqe3 := kq20;
kq20 = -7.19881671;
kqe2l := kq19;
kqe1l := kq18;
kqd4l := kq17;
kqd3l := kq16;
kqd2l := kq15;
kqd1l := kq14;
kqc3l := kq13;
kqc2l := kq12;
kqc1l := kq11;
kqb6l := kq10;
kqb5l := kq09;
kqb4l := kq08;
kqb3l := kq07;
kqb2l := kq06;
kqb1l := kq05;
kqa4l := kq04;
kqa3l := kq03;
kqa2l := kq02;
kqa1l := kq01;
mseqari: marker;
ibpm: monitor;
ibpma1: ibpm;
qa1r: quadrupole,l:= 0.21,k1:=kqa1r ;
qa2r: quadrupole,l:= 0.21,k1:=kqa2r ;
sqa1r: quadrupole,l:= 0.1,k1s:= 0;
ibpma2r: ibpm;
qa3r: quadrupole,l:= 0.21,k1:=kqa3r ;
qa4r: quadrupole,l:= 0.21,k1:=kqa4r ;
ibpma3r: ibpm;
sqa2r: quadrupole,l:= 0.1,k1s:= 0;
mseqare: marker;
mphm1ri: marker;
dedge30: dipedge,e1:= 0,h:= 1.338646717,hgap:= 0.015042,fint:= 0.5;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm1re: marker;
mseqbri: marker;
sqb1r: quadrupole,l:= 0.1,k1s:= 0;
qb1r: quadrupole,l:= 0.21,k1:=kqb1r ;
qb2r: quadrupole,l:= 0.21,k1:=kqb2r ;
qb3r: quadrupole,l:= 0.21,k1:=kqb3r ;
ibpmb1r: ibpm;
nlr1: marker;
ior: marker;
nlr2: marker;
ibpmb2r: ibpm;
qb4r: quadrupole,l:= 0.21,k1:=kqb4r ;
qb5r: quadrupole,l:= 0.21,k1:=kqb5r ;
qb6r: quadrupole,l:= 0.21,k1:=kqb6r ;
sqb2r: quadrupole,l:= 0.1,k1s:= 0;
mseqbre: marker;
mphm2ri: marker;
dedge60: dipedge,e1:= 0,h:= 1.381554029,hgap:= 0.014786,fint:= 0.5;
m2r: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm2re: marker;
mseqcri: marker;
sqc1r: quadrupole,l:= 0.1,k1s:= 0;
ibpmc1r: ibpm;
qc1r: quadrupole,l:= 0.21,k1:=kqc1r ;
qc2r: quadrupole,l:= 0.21,k1:=kqc2r ;
qc3r: quadrupole,l:= 0.21,k1:=kqc3r ;
ibpmc2r: ibpm;
sqc2r: quadrupole,l:= 0.1,k1s:= 0;
mseqcre: marker;
mphm3ri: marker;
m3r: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm3re: marker;
mseqdri: marker;
ibpmd1r: ibpm;
sqd1r: quadrupole,l:= 0.1,k1s:= 0;
qd1r: quadrupole,l:= 0.21,k1:=kqd1r ;
qd2r: quadrupole,l:= 0.21,k1:=kqd2r ;
el1: marker;
cel: solenoid,l:= 0.7,ks:= 0;
el2: marker;
qd3r: quadrupole,l:= 0.21,k1:=kqd3r ;
sqd2r: quadrupole,l:= 0.1,k1s:= 0;
qd4r: quadrupole,l:= 0.21,k1:=kqd4r ;
ibpmd2r: ibpm;
mseqdre: marker;
mphm4ri: marker;
m4r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm4re: marker;
mseqei: marker;
ibpme1r: ibpm;
qe1r: quadrupole,l:= 0.21,k1:=kqe1r ;
sqe1r: quadrupole,l:= 0.1,k1s:= 0;
ibpme2r: ibpm;
qe2r: quadrupole,l:= 0.21,k1:=kqe2r ;
sqe2r: quadrupole,l:= 0.1,k1s:= 0;
qe3: quadrupole,l:= 0.21,k1:=kqe3 ;
sqe2l: quadrupole,l:= 0.1,k1s:= 0;
qe2l: quadrupole,l:= 0.21,k1:=kqe2l ;
ibpme2l: ibpm;
sqe1l: quadrupole,l:= 0.1,k1s:= 0;
qe1l: quadrupole,l:= 0.21,k1:=kqe1l ;
ibpme1l: ibpm;
mseqee: marker;
mphm4li: marker;
m4l: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm4le: marker;
mseqdli: marker;
ibpmd2l: ibpm;
qd4l: quadrupole,l:= 0.21,k1:=kqd4l ;
sqd2l: quadrupole,l:= 0.1,k1s:= 0;
qd3l: quadrupole,l:= 0.21,k1:=kqd3l ;
rfc: rfcavity,l:= 0.05,volt:= 0.000847,lag:= 0,harmon:= 4;
qd2l: quadrupole,l:= 0.21,k1:=kqd2l ;
qd1l: quadrupole,l:= 0.21,k1:=kqd1l ;
sqd1l: quadrupole,l:= 0.1,k1s:= 0;
ibpmd1l: ibpm;
mseqdle: marker;
mphm3li: marker;
m3l: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm3le: marker;
mseqcli: marker;
sqc2l: quadrupole,l:= 0.1,k1s:= 0;
ibpmc2l: ibpm;
qc3l: quadrupole,l:= 0.21,k1:=kqc3l ;
qc2l: quadrupole,l:= 0.21,k1:=kqc2l ;
qc1l: quadrupole,l:= 0.21,k1:=kqc1l ;
ibpmc1l: ibpm;
sqc1l: quadrupole,l:= 0.1,k1s:= 0;
mseqcle: marker;
mphm2li: marker;
m2l: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm2le: marker;
mseqbli: marker;
sqb2l: quadrupole,l:= 0.1,k1s:= 0;
qb6l: quadrupole,l:= 0.21,k1:=kqb6l ;
qb5l: quadrupole,l:= 0.21,k1:=kqb5l ;
qb4l: quadrupole,l:= 0.21,k1:=kqb4l ;
ibpmb2l: ibpm;
nll1: marker;
nll2: marker;
ibpmb1l: ibpm;
qb3l: quadrupole,l:= 0.21,k1:=kqb3l ;
qb2l: quadrupole,l:= 0.21,k1:=kqb2l ;
qb1l: quadrupole,l:= 0.21,k1:=kqb1l ;
sqb1l: quadrupole,l:= 0.1,k1s:= 0;
mseqble: marker;
mphm1li: marker;
m1l: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm1le: marker;
mseqali: marker;
sqa2l: quadrupole,l:= 0.1,k1s:= 0;
ibpma3l: ibpm;
qa4l: quadrupole,l:= 0.21,k1:=kqa4l ;
qa3l: quadrupole,l:= 0.21,k1:=kqa3l ;
ibpma2l: ibpm;
sqa1l: quadrupole,l:= 0.1,k1s:= 0;
qa2l: quadrupole,l:= 0.21,k1:=kqa2l ;
qa1l: quadrupole,l:= 0.21,k1:=kqa1l ;
mseqale: marker;
iota: sequence, l = 39.95567226;
mseqari, at = 0;
ibpma1, at = 0.02;
qa1r, at = 1.0175;
qa2r, at = 1.3625;
sqa1r, at = 1.99;
ibpma2r, at = 2.095;
qa3r, at = 2.2975;
qa4r, at = 2.6525;
ibpma3r, at = 2.865;
sqa2r, at = 2.97;
mseqare, at = 3.0405;
mphm1ri, at = 3.0405;
dedge30, at = 3.117400202;
m1r, at = 3.312970389;
dedge30, at = 3.508540575;
mphm1re, at = 3.585440777;
mseqbri, at = 3.585440777;
sqb1r, at = 3.715940777;
qb1r, at = 3.953440777;
qb2r, at = 4.303440777;
qb3r, at = 4.653440777;
ibpmb1r, at = 4.865940777;
nlr1, at = 4.910940777;
ior, at = 5.810940777;
nlr2, at = 6.710940777;
ibpmb2r, at = 6.755940777;
qb4r, at = 6.968440777;
qb5r, at = 7.318440777;
qb6r, at = 7.668440777;
sqb2r, at = 7.905940777;
mseqbre, at = 8.036440777;
mphm2ri, at = 8.036440777;
dedge60, at = 8.112186805;
m2r, at = 8.491179421;
dedge60, at = 8.870172037;
mphm2re, at = 8.945918065;
mseqcri, at = 8.945918065;
sqc1r, at = 9.106418065;
ibpmc1r, at = 9.211418065;
qc1r, at = 9.423918065;
qc2r, at = 9.988918065;
qc3r, at = 10.55391806;
ibpmc2r, at = 10.76641806;
sqc2r, at = 10.87141806;
mseqcre, at = 11.03191806;
mphm3ri, at = 11.03191806;
dedge60, at = 11.10766409;
m3r, at = 11.48665671;
dedge60, at = 11.86564932;
mphm3re, at = 11.94139535;
mseqdri, at = 11.94139535;
ibpmd1r, at = 12.20689535;
sqd1r, at = 12.39189535;
qd1r, at = 12.61939535;
qd2r, at = 13.24939535;
el1, at = 13.81689535;
cel, at = 14.16689535;
el2, at = 14.51689535;
qd3r, at = 15.08439535;
sqd2r, at = 15.39939535;
qd4r, at = 15.71439535;
ibpmd2r, at = 16.12689535;
mseqdre, at = 16.39239535;
mphm4ri, at = 16.39239535;
dedge30, at = 16.46929555;
m4r, at = 16.66486574;
dedge30, at = 16.86043593;
mphm4re, at = 16.93733613;
mseqei, at = 16.93733613;
ibpme1r, at = 17.20283613;
qe1r, at = 17.51533613;
sqe1r, at = 17.76283613;
ibpme2r, at = 18.74783613;
qe2r, at = 18.96033613;
sqe2r, at = 19.34783613;
qe3, at = 19.97783613;
sqe2l, at = 20.60783613;
qe2l, at = 20.99533613;
ibpme2l, at = 21.20783613;
sqe1l, at = 22.19283613;
qe1l, at = 22.44033613;
ibpme1l, at = 22.75283613;
mseqee, at = 23.01833613;
mphm4li, at = 23.01833613;
dedge30, at = 23.09523633;
m4l, at = 23.29080652;
dedge30, at = 23.4863767;
mphm4le, at = 23.56327691;
mseqdli, at = 23.56327691;
ibpmd2l, at = 23.82877691;
qd4l, at = 24.24127691;
sqd2l, at = 24.55627691;
qd3l, at = 24.87127691;
rfc, at = 25.78877691;
qd2l, at = 26.70627691;
qd1l, at = 27.33627691;
sqd1l, at = 27.56377691;
ibpmd1l, at = 27.74877691;
mseqdle, at = 28.01427691;
mphm3li, at = 28.01427691;
dedge60, at = 28.09002293;
m3l, at = 28.46901555;
dedge60, at = 28.84800817;
mphm3le, at = 28.92375419;
mseqcli, at = 28.92375419;
sqc2l, at = 29.08425419;
ibpmc2l, at = 29.18925419;
qc3l, at = 29.40175419;
qc2l, at = 29.96675419;
qc1l, at = 30.53175419;
ibpmc1l, at = 30.74425419;
sqc1l, at = 30.84925419;
mseqcle, at = 31.00975419;
mphm2li, at = 31.00975419;
dedge60, at = 31.08550022;
m2l, at = 31.46449284;
dedge60, at = 31.84348545;
mphm2le, at = 31.91923148;
mseqbli, at = 31.91923148;
sqb2l, at = 32.04973148;
qb6l, at = 32.28723148;
qb5l, at = 32.63723148;
qb4l, at = 32.98723148;
ibpmb2l, at = 33.19973148;
nll1, at = 33.24473148;
nll2, at = 35.04473148;
ibpmb1l, at = 35.08973148;
qb3l, at = 35.30223148;
qb2l, at = 35.65223148;
qb1l, at = 36.00223148;
sqb1l, at = 36.23973148;
mseqble, at = 36.37023148;
mphm1li, at = 36.37023148;
dedge30, at = 36.44713168;
m1l, at = 36.64270187;
dedge30, at = 36.83827206;
mphm1le, at = 36.91517226;
mseqali, at = 36.91517226;
sqa2l, at = 36.98567226;
ibpma3l, at = 37.09067226;
qa4l, at = 37.30317226;
qa3l, at = 37.65817226;
ibpma2l, at = 37.86067226;
sqa1l, at = 37.96567226;
qa2l, at = 38.59317226;
qa1l, at = 38.93817226;
mseqale, at = 39.95567226;
endsequence;
beam, particle=proton, energy=pmass+0.0025;
"""


def test_linear(lattice, line):

    reader = synmadx.MadX_reader()
    reader.parse(lattice)

    lattice = reader.get_lattice(line)
    print("synergia lattice")
    print(lattice)
    print()

    ix_lattice = syn2_to_impactx(lattice, order=Order.linear)

    print("impactx lattice")
    print(ix_lattice)
    print()
    print("unrolled impactx lattice")
    print(unroll_impactx_lattice(ix_lattice))

    return


def test_exact(lattice, line):

    reader = synmadx.MadX_reader()
    reader.parse(lattice)

    lattice = reader.get_lattice(line)
    print("synergia lattice")
    print(lattice)
    print()

    ix_lattice = syn2_to_impactx(lattice, order=Order.exact)

    print("impactx lattice")
    print(ix_lattice)
    print()
    print("unrolled impactx lattice")
    print(unroll_impactx_lattice(ix_lattice))

    return


if __name__ == "__main__":
    # test_linear(__misc_txt, "misc")
    test_exact(__misc_txt, "misc")

    # test_linear(__simple_booster_txt, "booster")
    # test_exact(__simple_booster_txt, "booster")

    # test_linear(__iota_txt, "iota")
    # test_exact(__iota_txt, "iota")
