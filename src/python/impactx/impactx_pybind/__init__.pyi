"""

impactx_pybind
--------------
.. currentmodule:: impactx_pybind

.. autosummary::
   :toctree: _generate
   ImpactX
   distribution
   elements

"""

from __future__ import annotations

import typing

import pybind11_stubgen.typing_ext

import amrex.space3d.amrex_3d_pybind
from amrex import space3d as amr

from . import distribution, elements, wakeconvolution

__all__ = [
    "Config",
    "CoordSystem",
    "ImpactX",
    "ImpactXParConstIter",
    "ImpactXParIter",
    "ImpactXParticleContainer",
    "RefPart",
    "amr",
    "coordinate_transformation",
    "distribution",
    "elements",
    "push",
    "s",
    "t",
    "wakeconvolution",
]

class Config:
    gpu_backend = None
    have_gpu: typing.ClassVar[bool] = False
    have_mpi: typing.ClassVar[bool] = True
    have_omp: typing.ClassVar[bool] = True

class CoordSystem:
    """
    Members:

      s

      t
    """

    __members__: typing.ClassVar[
        dict[str, CoordSystem]
    ]  # value = {'s': <CoordSystem.s: 0>, 't': <CoordSystem.t: 1>}
    s: typing.ClassVar[CoordSystem]  # value = <CoordSystem.s: 0>
    t: typing.ClassVar[CoordSystem]  # value = <CoordSystem.t: 1>
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class ImpactX:
    def DistributionMap(
        self, lev: int
    ) -> amrex.space3d.amrex_3d_pybind.DistributionMapping: ...
    def Geom(self, lev: int) -> amrex.space3d.amrex_3d_pybind.Geometry: ...
    def __init__(self) -> None: ...
    def add_particles(
        self,
        bunch_charge: float,
        distr: distribution.Empty
        | distribution.Gaussian
        | distribution.Kurth4D
        | distribution.Kurth6D
        | distribution.KVdist
        | distribution.Thermal
        | distribution.Triangle
        | distribution.Semigaussian
        | distribution.Waterbag,
        npart: int,
    ) -> None:
        """
        Generate and add n particles to the particle container.

        Will also resize the geometry based on the updated particle
        distribution's extent and then redistribute particles in according
        AMReX grid boxes.
        """
    def boxArray(self, lev: int) -> amrex.space3d.amrex_3d_pybind.BoxArray: ...
    def deposit_charge(self) -> None:
        """
        Deposit charge in x,y,z.
        """
    def evolve(self) -> None:
        """
        Run the main simulation loop for a number of steps.
        """
    def finalize(self) -> None:
        """
        Deallocate all contexts and data.
        """
    def init_beam_distribution_from_inputs(self) -> None: ...
    def init_grids(self) -> None:
        """
        Initialize AMReX blocks/grids for domain decomposition & space charge mesh.

        This must come first, before particle beams and lattice elements are initialized.
        """
    def init_lattice_elements_from_inputs(self) -> None: ...
    def load_inputs_file(self, arg0: str) -> None: ...
    def particle_container(self) -> ImpactXParticleContainer:
        """
        Access the beam particle container.
        """
    def phi(self, lev: int) -> amrex.space3d.amrex_3d_pybind.MultiFab:
        """
        scalar potential per level
        """
    def resize_mesh(self) -> None:
        """
        Resize the mesh :py:attr:`~domain` based on the :py:attr:`~dynamic_size` and related parameters.
        """
    def rho(self, lev: int) -> amrex.space3d.amrex_3d_pybind.MultiFab:
        """
        charge density per level
        """
    def space_charge_field(
        self, lev: int, comp: str
    ) -> amrex.space3d.amrex_3d_pybind.MultiFab:
        """
        space charge force (vector: x,y,z) per level
        """
    @property
    def abort_on_unused_inputs(self) -> int:
        """
        Configure simulation to abort AFTER it has run
        if there are unused parameters in the input.
        """
    @abort_on_unused_inputs.setter
    def abort_on_unused_inputs(self, arg1: int) -> None: ...
    @property
    def abort_on_warning_threshold(self) -> str:
        """
        Set WarnPriority threshold to decide if ImpactX
        has to abort when a warning is recorded.
        Valid choices are: ['low', 'medium', 'high'].
        """
    @abort_on_warning_threshold.setter
    def abort_on_warning_threshold(self, arg1: str) -> None: ...
    @property
    def always_warn_immediately(self) -> int:
        """
        If set to 1, immediately prints every warning message
         as soon as it is generated.
        """
    @always_warn_immediately.setter
    def always_warn_immediately(self, arg1: int) -> None: ...
    @property
    def blocking_factor_x(self) -> list[int]:
        """
        AMReX blocking factor for a direction, per MR level.
        """
    @blocking_factor_x.setter
    def blocking_factor_x(self, arg1: list[int]) -> None: ...
    @property
    def blocking_factor_y(self) -> list[int]:
        """
        AMReX blocking factor for a direction, per MR level.
        """
    @blocking_factor_y.setter
    def blocking_factor_y(self, arg1: list[int]) -> None: ...
    @property
    def blocking_factor_z(self) -> list[int]:
        """
        AMReX blocking factor for a direction, per MR level.
        """
    @blocking_factor_z.setter
    def blocking_factor_z(self, arg1: list[int]) -> None: ...
    @property
    def csr(self) -> bool:
        """
        Enable or disable Coherent Synchrotron Radiation (CSR) calculations (default: disabled).
        """
    @csr.setter
    def csr(self, arg1: bool) -> None: ...
    @property
    def csr_bins(self) -> bool:
        """
        Number of longitudinal bins used for CSR calculations (default: 150).
        """
    @csr_bins.setter
    def csr_bins(self, arg1: int) -> None: ...
    @property
    def diag_file_min_digits(self) -> int:
        """
        The minimum number of digits (default: 6) used for the step
        number appended to the diagnostic file names.
        """
    @diag_file_min_digits.setter
    def diag_file_min_digits(self, arg1: int) -> None: ...
    @property
    def diagnostics(self) -> bool:
        """
        Enable or disable diagnostics generally (default: enabled).
        Disabling this is mostly used for benchmarking.
        """
    @diagnostics.setter
    def diagnostics(self, arg1: bool) -> None: ...
    @property
    def domain(self) -> amrex.space3d.amrex_3d_pybind.RealBox:
        """
        The physical extent of the full simulation domain, relative to the reference particle position, in meters.
        """
    @domain.setter
    def domain(self, arg1: amrex.space3d.amrex_3d_pybind.RealBox) -> None: ...
    @property
    def dynamic_size(self) -> bool:
        """
        Use dynamic (``true``) resizing of the field mesh or static sizing (``false``).
        """
    @dynamic_size.setter
    def dynamic_size(self, arg1: bool) -> None: ...
    @property
    def finest_level(self) -> int:
        """
        The currently finest level of mesh-refinement used. This is always less or equal to max_level.
        """
    @property
    def lattice(self) -> elements.KnownElementsList:
        """
        Access the accelerator element lattice.
        """
    @lattice.setter
    def lattice(self, arg0: elements.KnownElementsList) -> None: ...
    @property
    def max_level(self) -> int:
        """
        The maximum mesh-refinement level for the simulation.
        """
    @max_level.setter
    def max_level(self, arg1: int) -> None: ...
    @property
    def mlmg_absolute_tolerance(self) -> bool:
        """
        The absolute tolerance with which the space-charge fields should be calculated in units of V/m^2. More specifically, the acceptable residual with which the solution can be considered converged. In general this should be left as the default, but in cases where the simulation state changes very little between steps it can occur that the initial guess for the MLMG solver is so close to the converged value that it fails to improve that solution sufficiently to reach the mlmg_relative_tolerance value.
        """
    @mlmg_absolute_tolerance.setter
    def mlmg_absolute_tolerance(self, arg1: float) -> None: ...
    @property
    def mlmg_max_iters(self) -> bool:
        """
        Maximum number of iterations used for MLMG solver for space-charge fields calculation. In case if MLMG converges but fails to reach the desired self_fields_required_precision, this parameter may be increased.
        """
    @mlmg_max_iters.setter
    def mlmg_max_iters(self, arg1: int) -> None: ...
    @property
    def mlmg_relative_tolerance(self) -> bool:
        """
        The relative precision with which the electrostatic space-charge fields should be calculated. More specifically, the space-charge fields are computed with an iterative Multi-Level Multi-Grid (MLMG) solver. This solver can fail to reach the default precision within a reasonable time.
        """
    @mlmg_relative_tolerance.setter
    def mlmg_relative_tolerance(self, arg1: float) -> None: ...
    @property
    def mlmg_verbosity(self) -> bool:
        """
        The verbosity used for MLMG solver for space-charge fields calculation. Currently MLMG solver looks for verbosity levels from 0-5. A higher number results in more verbose output.
        """
    @mlmg_verbosity.setter
    def mlmg_verbosity(self, arg1: int) -> None: ...
    @property
    def n_cell(self) -> list[int]:
        """
        The number of grid points along each direction on the coarsest level.
        """
    @n_cell.setter
    def n_cell(
        self,
        arg1: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)],
    ) -> None: ...
    @property
    def particle_lost_diagnostics_backend(self) -> str:
        """
        Diagnostics for particles lost in apertures.

        See the ``BeamMonitor`` element for backend values.
        """
    @particle_lost_diagnostics_backend.setter
    def particle_lost_diagnostics_backend(self, arg1: str) -> None: ...
    @property
    def particle_shape(self) -> int:
        """
        Whether to calculate space charge effects.
        """
    @particle_shape.setter
    def particle_shape(self, arg1: int) -> None: ...
    @property
    def periods(self) -> int:
        """
        The number of periods to repeat the lattice.
        """
    @periods.setter
    def periods(self, arg1: int) -> None: ...
    @property
    def poisson_solver(self) -> str:
        """
        The numerical solver to solve the Poisson equation when calculating space charge effects. Either multigrid (default) or fft.
        """
    @poisson_solver.setter
    def poisson_solver(self, arg1: str) -> None: ...
    @property
    def prob_relative(self) -> float:
        """
        The field mesh spans, per direction, multiple times the maximum physical extent of beam particles, as given by this factor.
        """
    @prob_relative.setter
    def prob_relative(self, arg1: list[float]) -> None: ...
    @property
    def slice_step_diagnostics(self) -> bool:
        """
        Enable or disable diagnostics every slice step in elements (default: disabled).

        By default, diagnostics is performed at the beginning and end of the simulation.
        Enabling this flag will write diagnostics every step and slice step.
        """
    @slice_step_diagnostics.setter
    def slice_step_diagnostics(self, arg1: bool) -> None: ...
    @property
    def space_charge(self) -> bool:
        """
        Enable or disable space charge calculations (default: enabled).
        """
    @space_charge.setter
    def space_charge(self, arg1: bool) -> None: ...
    @property
    def verbose(self) -> int:
        """
        Controls how much information is printed to the terminal, when running ImpactX.
        ``0`` for silent, higher is more verbose. Default is ``1``.
        """
    @verbose.setter
    def verbose(self, arg1: int) -> None: ...

class ImpactXParConstIter(
    amrex.space3d.amrex_3d_pybind.ParConstIter_pureSoA_8_0_default
):
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_8_0_default,
        level: int,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_8_0_default,
        level: int,
        info: amrex.space3d.amrex_3d_pybind.MFItInfo,
    ) -> None: ...
    def pc(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_8_0_default: ...
    def soa(self):
        """
        Get the StructOfArrays on the current tile

            Parameters
            ----------
            self : ImpactXParIter or ImpactXParConstIter
              used to query particle container component names

        """

class ImpactXParIter(amrex.space3d.amrex_3d_pybind.ParIter_pureSoA_8_0_default):
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_8_0_default,
        level: int,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_8_0_default,
        level: int,
        info: amrex.space3d.amrex_3d_pybind.MFItInfo,
    ) -> None: ...
    def pc(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_8_0_default: ...
    def soa(self):
        """
        Get the StructOfArrays on the current tile

            Parameters
            ----------
            self : ImpactXParIter or ImpactXParConstIter
              used to query particle container component names

        """

class ImpactXParticleContainer(
    amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_8_0_default
):
    const_iterator = ImpactXParConstIter
    iterator = ImpactXParIter
    def add_n_particles(
        self,
        x: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        y: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        t: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        px: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        py: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        pt: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        qm: float,
        bchchg: float,
    ) -> None:
        """
        Add new particles to the container for fixed s.

        Note: This can only be used *after* the initialization (grids) have
              been created, meaning after the call to ImpactX.init_grids
              has been made in the ImpactX class.

        :param x: positions in x
        :param y: positions in y
        :param t: positions as time-of-flight in c*t
        :param px: momentum in x
        :param py: momentum in y
        :param pt: momentum in t
        :param qm: charge over mass in 1/eV
        :param bchchg: total charge within a bunch in C
        """
    def mean_and_std_positions(self) -> tuple[float, float, float, float, float, float]:
        """
        Compute the mean and std of the particle position in each dimension.

        :return: x_mean, x_std, y_mean, y_std, z_mean, z_std
        """
    def min_and_max_positions(self) -> tuple[float, float, float, float, float, float]:
        """
        Compute the min and max of the particle position in each dimension.

        :return: x_min, y_min, z_min, x_max, y_max, z_max
        """
    def plot_phasespace(self, num_bins=50, root_rank=0):
        """

        Plot the longitudinal and transverse phase space projections with matplotlib.

        Parameters
        ----------
        self : ImpactXParticleContainer_*
            The particle container class in ImpactX
        num_bins : int, default=50
            The number of bins for spatial and momentum directions per plot axis.
        root_rank : int, default=0
            MPI root rank to reduce to in parallel runs.

        Returns
        -------
        A matplotlib figure with containing the plot.
        For MPI-parallel ranks, the figure is only created on the root_rank.

        """
    def redistribute(
        self, arg0: int, arg1: int, arg2: int, arg3: int, arg4: bool
    ) -> None:
        """
        Redistribute particles in the current mesh in x, y, z
        """
    def reduced_beam_characteristics(self) -> dict[str, float]:
        """
        Compute reduced beam characteristics like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters.
        """
    def ref_particle(self) -> RefPart:
        """
        Access the reference particle.
        """
    def set_ref_particle(self, refpart: RefPart) -> None:
        """
        Set reference particle attributes.
        """
    @property
    def RealSoA_names(self) -> list[str]:
        """
        Get the name of each ParticleReal SoA component
        """
    @property
    def coord_system(self) -> CoordSystem:
        """
        Get the current coordinate system of particles in this container
        """
    @property
    def intSoA_names(self) -> list[str]:
        """
        Get the name of each int SoA component
        """

class RefPart:
    @staticmethod
    def load_file(ref: RefPart, madx_file):
        """

        Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
        :param RefPart ref: ImpactX reference particle (passed by reference)
        :param madx_file: file name to MAD-X file with beamline elements
        :return: list of ImpactX.KnownElements

        """
    def __init__(self) -> None:
        """
        This struct stores the reference particle attributes
        stored in ImpactXParticleContainer.
        """
    def set_charge_qe(self, charge_qe: float) -> RefPart:
        """
        Set reference particle charge (positive elementary charge)
        """
    def set_kin_energy_MeV(self, kin_energy_MeV: float) -> RefPart:
        """
        Set reference particle kinetic energy (MeV)
        """
    def set_mass_MeV(self, mass_MeV: float) -> RefPart:
        """
        Set reference particle rest mass (MeV/c^2)
        """
    @property
    def beta(self) -> float:
        """
        Get reference particle relativistic beta
        """
    @property
    def beta_gamma(self) -> float:
        """
        Get reference particle beta*gamma
        """
    @property
    def charge(self) -> float:
        """
        reference charge, in C
        """
    @charge.setter
    def charge(self, arg0: float) -> None: ...
    @property
    def charge_qe(self) -> float:
        """
        Get reference particle charge (positive elementary charge)
        """
    @property
    def gamma(self) -> float:
        """
        Get reference particle relativistic gamma
        """
    @property
    def kin_energy_MeV(self) -> float:
        """
        Get reference particle energy (MeV)
        """
    @property
    def mass(self) -> float:
        """
        reference rest mass, in kg
        """
    @mass.setter
    def mass(self, arg0: float) -> None: ...
    @property
    def mass_MeV(self) -> float:
        """
        Get reference particle rest mass (MeV/c^2)
        """
    @property
    def pt(self) -> float:
        """
        energy deviation, normalized by rest energy
        """
    @pt.setter
    def pt(self, arg0: float) -> None: ...
    @property
    def px(self) -> float:
        """
        momentum in x, normalized to proper velocity
        """
    @px.setter
    def px(self, arg0: float) -> None: ...
    @property
    def py(self) -> float:
        """
        momentum in y, normalized to proper velocity
        """
    @py.setter
    def py(self, arg0: float) -> None: ...
    @property
    def pz(self) -> float:
        """
        momentum in z, normalized to proper velocity
        """
    @pz.setter
    def pz(self, arg0: float) -> None: ...
    @property
    def qm_ratio_SI(self) -> float:
        """
        Get reference particle charge to mass ratio (C/kg)
        """
    @property
    def rigidity_Tm(self) -> float:
        """
        Get reference particle magnetic rigidity Brho (T*m)
        """
    @property
    def s(self) -> float:
        """
        integrated orbit path length, in meters
        """
    @s.setter
    def s(self, arg0: float) -> None: ...
    @property
    def t(self) -> float:
        """
        clock time * c in meters
        """
    @t.setter
    def t(self, arg0: float) -> None: ...
    @property
    def x(self) -> float:
        """
        horizontal position x, in meters
        """
    @x.setter
    def x(self, arg0: float) -> None: ...
    @property
    def y(self) -> float:
        """
        vertical position y, in meters
        """
    @y.setter
    def y(self, arg0: float) -> None: ...
    @property
    def z(self) -> float:
        """
        longitudinal position y, in meters
        """
    @z.setter
    def z(self, arg0: float) -> None: ...

def coordinate_transformation(
    pc: ImpactXParticleContainer, direction: CoordSystem
) -> None:
    """
    Transform coordinates from fixed s to fixed to or vice versa.
    """

def push(
    pc: ImpactXParticleContainer,
    element: elements.Empty
    | elements.Aperture
    | elements.Buncher
    | elements.CFbend
    | elements.ChrAcc
    | elements.ChrDrift
    | elements.ChrPlasmaLens
    | elements.ChrQuad
    | elements.ConstF
    | elements.BeamMonitor
    | elements.DipEdge
    | elements.Drift
    | elements.ExactDrift
    | elements.ExactSbend
    | elements.Kicker
    | elements.Multipole
    | elements.NonlinearLens
    | elements.Programmable
    | elements.PRot
    | elements.Quad
    | elements.RFCavity
    | elements.Sbend
    | elements.ShortRF
    | elements.SoftSolenoid
    | elements.SoftQuadrupole
    | elements.Sol
    | elements.TaperedPL
    | elements.ThinDipole,
    step: int = 0,
) -> None:
    """
    Push particles through an element
    """

__author__: str = (
    "Axel Huebl, Chad Mitchell, Ryan Sandberg, Marco Garten, Ji Qiang, et al."
)
__license__: str = "BSD-3-Clause-LBNL"
__version__: str = "24.08"
s: CoordSystem  # value = <CoordSystem.s: 0>
t: CoordSystem  # value = <CoordSystem.t: 1>
