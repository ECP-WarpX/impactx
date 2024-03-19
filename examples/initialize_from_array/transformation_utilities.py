import numpy as np

def to_ref_part_t_from_lab_t(ref_part, x, y, z, px, py, pz):
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
def from_ref_part_t_to_lab_t(ref_part, dx, dy, dz, dpx, dpy, dpz):
    x = dx + ref_part.x
    y = dy + ref_part.y
    z = dz + ref_part.z

    px = ref_part.px + ref_part.pz * dpx
    py = ref_part.py + ref_part.pz * dpy
    pz = ref_part.pz + ref_part.pz * dpz

    return x, y, z, px, py, pz
def to_s_from_t(ref_part, dx,dy,dz,dpx,dpy,dpz): #data_arr_t):
    """
    """
    ref_pz = ref_part.pz
    ref_pt = ref_part.pt
    dxs = dx - ref_pz*dpx*dz/(ref_pz+ref_pz*dpz)
    dys = dy - ref_pz*dpy*dz/(ref_pz+ref_pz*dpz)
    pt = -np.sqrt(1 + (ref_pz + ref_pz*dpz)**2+(ref_pz*dpx)**2+(ref_pz*dpy)**2)
    dt = pt*dz/(ref_pz + ref_pz*dpz)
    dpt = (pt - ref_pt) / ref_pz
    return dxs, dys, dt, dpx, dpy, dpt
def to_t_from_s(ref_part, dx,dy,dt,dpx,dpy,dpt): #data_arr_t):
    """
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
    def __init__(self,x,y,z,px,py,pz, pt):
        self.attr_list = ['x','y','z','px','py','pz','pt']
        self.x = x
        self.y = y
        self.z = z
        self.px = px
        self.py = py
        self.pz = pz
        self.pt = pt

    def __repr__(self):
        mystr = ''
        for attr in self.attr_list:
            mystr += f'self.{attr}={getattr(self,attr)}, ' 
        return mystr
    

    def __str__(self):
        mystr = ''
        for attr in self.attr_list:
            mystr += f'self.{attr}={getattr(self,attr)}, '
        return mystr

    