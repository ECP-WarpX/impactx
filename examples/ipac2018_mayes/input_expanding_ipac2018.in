###############################################################################
# Particle Beam(s)
###############################################################################
beam.npart = 10000  # outside tests, use 1e5 or more
beam.units = static
beam.energy = 9.48900105  # KE corresponding to 10 MeV total energy
beam.charge = 1.0e-9  # 1 nC
beam.particle = electron
beam.distribution = gaussian
beam.sigmaX = 1.0e-3  # 1 mm
beam.sigmaY = 1.0e-3  # 1 mm
beam.sigmaT = 1.00130816210e-4  # corresponding to sig_z = 0.1 mm
beam.sigmaPx = 0.0
beam.sigmaPy = 0.0
beam.sigmaPt = 0.0


###############################################################################
# Beamline: lattice elements and segments
###############################################################################
lattice.elements = monitor drift1 monitor
lattice.nslice = 40

drift1.type = drift
drift1.ds = 1.0  # 1 m drift

monitor.type = beam_monitor
monitor.backend = h5


###############################################################################
# Algorithms
###############################################################################
algo.particle_shape = 2
algo.space_charge = true

amr.n_cell = 56 56 48
geometry.prob_relative = 3.0
