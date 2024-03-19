import numpy as np


def to_ref_part_t_from_global_t(ref_part, x, y, z, px, py, pz):
    """
    Transform from global coordinates to reference particle coordinates
    
    This function takes particle bunch data and returns the bunch phase space positions relative to a reference particle

    Parameters
    ---
    ref_part: reference particle object, either from ImpactX or from the class defined in this file, MyRefPart
    x: array-like, global beam particle x-positions
    y: array-like, global beam particle y-positions
    z: array-like, global beam particle z-positions
    px: array-like, global beam particle x-momenta
    py: array-like, global beam particle y-momenta
    pz: array-like, global beam particle z-momenta

    Returns
    -------
    dx: array-like, beam particle x-positions relative to reference particle x value, ``ref_part.x``
    dy: array-like, beam particle y-positions relative to reference particle y value, ``ref_part.y``
    dz: array-like, beam particle z-positions relative to reference particle z value, ``ref_part.z``
    dpx: array-like, beam particle x-momenta relative to reference particle x value, ``ref_part.px``
    dpy: array-like, beam particle y-momenta relative to reference particle y value, ``ref_part.py``
    dpz: array-like, beam particle z-momenta relative to reference particle z value, ``ref_part.pz``
    """
    dx = x - ref_part.x
    dy = y - ref_part.y
    dz = z - ref_part.z
    dpx = px - ref_part.px
    dpy = py - ref_part.py
    dpz = pz - ref_part.pz

    dpx /= ref_part.pz
    dpy /= ref_part.pz
    dpz /= ref_part.pz

    return dx, dy, dz, dpx, dpy, dpz


def to_global_t_from_ref_part_t(ref_part, dx, dy, dz, dpx, dpy, dpz):
    """
    Transform from reference particle to global coordinates
    
    This function takes particle bunch data relative to a reference particle
    and returns all particle data in the global coordinate frame

    Parameters
    ---
    ref_part: reference particle object, either from ImpactX or from the class defined in this file, MyRefPart
    dx: array-like, beam particle x-positions relative to reference particle x value, ``ref_part.x``
    dy: array-like, beam particle y-positions relative to reference particle y value, ``ref_part.y``
    dz: array-like, beam particle z-positions relative to reference particle z value, ``ref_part.z``
    dpx: array-like, beam particle x-momenta relative to reference particle x value, ``ref_part.px``
    dpy: array-like, beam particle y-momenta relative to reference particle y value, ``ref_part.py``
    dpz: array-like, beam particle z-momenta relative to reference particle z value, ``ref_part.pz``

    Returns
    -------
    x: array-like, global beam particle x-positions
    y: array-like, global beam particle y-positions
    z: array-like, global beam particle z-positions
    px: array-like, global beam particle x-momenta
    py: array-like, global beam particle y-momenta
    pz: array-like, global beam particle z-momenta
    """
    x = dx + ref_part.x
    y = dy + ref_part.y
    z = dz + ref_part.z

    px = ref_part.px + ref_part.pz * dpx
    py = ref_part.py + ref_part.pz * dpy
    pz = ref_part.pz + ref_part.pz * dpz

    return x, y, z, px, py, pz


def to_s_from_t(ref_part, dx, dy, dz, dpx, dpy, dpz):  # data_arr_t):
    """
    Transform from fixed-s to fixed-t coordinates

    This function takes particle bunch data relative to a reference particle
    at a fixed time t and returns data at a fixed longitudinal position s.
    That is, spatial distance is transformed to time delay from the reference particle.

    Parameters
    ---
    ref_part: reference particle object, either from ImpactX or from the class defined in this file, MyRefPart
    dx: array-like, beam particle x-positions relative to reference particle at fixed t
    dy: array-like, beam particle y-positions relative to reference particle at fixed t
    dz: array-like, beam particle z-positions relative to reference particle at fixed t
    dpx: array-like, beam particle x-momenta relative to reference particle at fixed t
    dpy: array-like, beam particle y-momenta relative to reference particle at fixed t
    dpz: array-like, beam particle z-momenta relative to reference particle at fixed t

    Returns
    -------
    dxs: array-like, beam particle x-positions relative to reference particle at fixed s
    dys: array-like, beam particle y-positions relative to reference particle at fixed s
    dt: array-like, beam particle time delay relative to reference particle at fixed s
    dpx: array-like, beam particle x-momenta relative to reference particle at fixed s
    dpy: array-like, beam particle y-momenta relative to reference particle at fixed s
    dpz: array-like, beam particle t-momenta (-gamma) relative to reference particle at fixed s
    """
    ref_pz = ref_part.pz
    ref_pt = ref_part.pt
    dxs = dx - ref_pz * dpx * dz / (ref_pz + ref_pz * dpz)
    dys = dy - ref_pz * dpy * dz / (ref_pz + ref_pz * dpz)
    pt = -np.sqrt(
        1 + (ref_pz + ref_pz * dpz) ** 2 + (ref_pz * dpx) ** 2 + (ref_pz * dpy) ** 2
    )
    dt = pt * dz / (ref_pz + ref_pz * dpz)
    dpt = (pt - ref_pt) / ref_pz
    return dxs, dys, dt, dpx, dpy, dpt


def to_t_from_s(ref_part, dx, dy, dt, dpx, dpy, dpt):  # data_arr_t):
    """ 
    Transform from fixed-t to fixed-s coordinates

    This function takes particle bunch data relative to a reference particle
    at a fixed longitudinal position s and returns data at a fixed time t.
    That is, time delay is transformed to spatial distance from the reference particle.

    Parameters
    ---
    ref_part: reference particle object, either from ImpactX or from the class defined in this file, MyRefPart
    dx: array-like, beam particle x-positions relative to reference particle at fixed s
    dy: array-like, beam particle y-positions relative to reference particle at fixed s
    dz: array-like, beam particle time delay relative to reference particle at fixed s
    dpx: array-like, beam particle x-momenta relative to reference particle at fixed s
    dpy: array-like, beam particle y-momenta relative to reference particle at fixed s
    dpt: array-like, beam particle t-momenta (-gamma) relative to reference particle at fixed s

    Returns
    -------
    dxt: array-like, beam particle x-positions relative to reference particle at fixed t
    dyt: array-like, beam particle y-positions relative to reference particle at fixed t
    dt: array-like, beam particle z-positions relative to reference particle at fixed t
    dpx: array-like, beam particle x-momenta relative to reference particle at fixed t
    dpy: array-like, beam particle y-momenta relative to reference particle at fixed t
    dpz: array-like, beam particle z-momenta relative to reference particle at fixed t
    """
    ref_pz = ref_part.pz
    ref_pt = ref_part.pt
    dxt = dx + ref_pz * dpx * dt / (ref_pt + ref_pz * dpt)
    dyt = dy + ref_pz * dpy * dt / (ref_pt + ref_pz * dpt)

    pz = np.sqrt(
        -1 + (ref_pt + ref_pz * dpt) ** 2 - (ref_pz * dpx) ** 2 - (ref_pz * dpy) ** 2
    )
    dz = dt * pz / (ref_pt + ref_pz * dpt)
    dpz = (pz - ref_pz) / ref_pz
    return dxt, dyt, dz, dpx, dpy, dpz


class MyRefPart:
    """Struct containing reference particle data
    
    This class replicates the data storage of the ImpactX reference particle.
    It is used in coordinate transformations when an ImpactX reference particle isn't available,
    so the transformation syntax works in the context of pyImpactX
    """
    def __init__(self, x, y, z, px, py, pz, pt):
        self.attr_list = ["x", "y", "z", "px", "py", "pz", "pt"]
        self.x = x
        self.y = y
        self.z = z
        self.px = px
        self.py = py
        self.pz = pz
        self.pt = pt

    def __repr__(self):
        mystr = ""
        for attr in self.attr_list:
            mystr += f"self.{attr}={getattr(self,attr)}, "
        return mystr

    def __str__(self):
        mystr = ""
        for attr in self.attr_list:
            mystr += f"self.{attr}={getattr(self,attr)}, "
        return mystr
