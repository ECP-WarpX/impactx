import numpy as np

try:
    import torch

    on_gpu = torch.cuda.is_available()
except ImportError:
    on_gpu = False

import amrex.space3d as amr
import transformation_utilities as pycoord
from impactx import ImpactX, elements

################

N_part = int(1e5)
beam_radius = 2e-3
sigr = 500e-6
sigpx = 10
sigpy = 10
px = np.random.normal(0, sigpx, N_part)
py = np.random.normal(0, sigpy, N_part)
theta = 2 * np.pi * np.random.rand(N_part)
r = abs(np.random.normal(beam_radius, sigr, N_part))
x = r * np.cos(theta)
y = r * np.sin(theta)
z_mean = 0
pz_mean = 2e4
z_std = 1e-3
pz_std = 2e2
zpz_std = -0.18
zpz_cov_list = [[z_std**2, zpz_std], [zpz_std, pz_std**2]]
z, pz = np.random.default_rng().multivariate_normal([0, 0], zpz_cov_list, N_part).T
pz += pz_mean

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

energy_gamma = np.sqrt(1 + pz_mean**2)
energy_MeV = 0.510998950 * energy_gamma  # reference energy
bunch_charge_C = 10.0e-15  # used with space charge

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(energy_MeV)
qm_eev = -1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV
ref.z = 0

pc = sim.particle_container()

dx, dy, dz, dpx, dpy, dpz = pycoord.to_ref_part_t_from_global_t(ref, x, y, z, px, py, pz)
dx, dy, dt, dpx, dpy, dpt = pycoord.to_s_from_t(ref, dx, dy, dz, dpx, dpy, dpz)

if not on_gpu:  # initialize using cpu-based PODVectors
    dx_podv = amr.PODVector_real_std()
    dy_podv = amr.PODVector_real_std()
    dt_podv = amr.PODVector_real_std()
    dpx_podv = amr.PODVector_real_std()
    dpy_podv = amr.PODVector_real_std()
    dpt_podv = amr.PODVector_real_std()
else:  # initialize on device using arena/gpu-based PODVectors
    dx_podv = amr.PODVector_real_arena()
    dy_podv = amr.PODVector_real_arena()
    dt_podv = amr.PODVector_real_arena()
    dpx_podv = amr.PODVector_real_arena()
    dpy_podv = amr.PODVector_real_arena()
    dpt_podv = amr.PODVector_real_arena()

for element in dx:
    dx_podv.push_back(element)
for element in dy:
    dy_podv.push_back(element)
for element in dt:
    dt_podv.push_back(element)
for element in dpx:
    dpx_podv.push_back(element)
for element in dpy:
    dpy_podv.push_back(element)
for element in dpt:
    dpt_podv.push_back(element)

pc.add_n_particles(
    dx_podv, dy_podv, dt_podv, dpx_podv, dpy_podv, dpt_podv, qm_eev, bunch_charge_C
)

monitor = elements.BeamMonitor("monitor")
sim.lattice.extend(
    [
        monitor,
        elements.Drift(ds=0.01),
        monitor,
    ]
)

sim.evolve()

# clean shutdown
sim.finalize()
