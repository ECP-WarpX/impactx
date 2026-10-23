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

import collections.abc
import typing

import amrex.space3d.amrex_3d_pybind
from amrex import space3d as amr
from amrex.space3d.amrex_3d_pybind import SmallMatrix_3x1_F_SI1_double as Vector3
from amrex.space3d.amrex_3d_pybind import SmallMatrix_3x6_F_SI1_double as Map3x6
from amrex.space3d.amrex_3d_pybind import SmallMatrix_6x6_F_SI1_double as Map6x6

from . import distribution, elements, wakeconvolution

__all__: list[str] = [
    "Config",
    "CoordSystem",
    "Envelope",
    "ImpactX",
    "ImpactXParConstIter",
    "ImpactXParIter",
    "ImpactXParticleContainer",
    "Map3x6",
    "Map6x6",
    "RefPart",
    "UnorderedMap",
    "Vector3",
    "amr",
    "coordinate_transformation",
    "create_envelope",
    "distribution",
    "elements",
    "flatten_charge_to_2D",
    "push",
    "reverse",
    "s",
    "t",
    "wakeconvolution",
]

class RefPart:
    @staticmethod
    def load_file(ref: RefPart, madx_file):
        """
        Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements

        .. warning::

           Our MAD-X parser is under active development and provided
           as a preview. Please check any loaded MAD-X beams very
           carefully. Please report your experience and bugs on our
           `issue tracker <https://github.com/BLAST-ImpactX/impactx/issues>`__.

        :param RefPart ref: ImpactX reference particle (passed by reference)
        :param madx_file: file name to MAD-X file with beamline elements
        :return: list of ImpactX.KnownElements
        """
    def __init__(self) -> None:
        """
        This struct stores the reference particle attributes
        stored in ImpactXParticleContainer.
        """
    def copy(self) -> RefPart:
        """
        Copy the reference particle
        """
    def reset(self, keep_mass: bool = False, keep_charge: bool = False) -> None:
        """
        Reset the reference particle
        """
    def set_charge_qe(
        self, charge_qe: typing.SupportsFloat | typing.SupportsIndex
    ) -> RefPart:
        """
        Set reference particle charge (positive elementary charge) [q_e]
        """
    def set_gyromagnetic_anomaly(
        self, gyromagnetic_anomaly: typing.SupportsFloat | typing.SupportsIndex
    ) -> RefPart:
        """
        Set reference particle gyromagnetic anomaly value (for spin tracking)
        """
    def set_kin_energy_MeV(
        self, kin_energy_MeV: typing.SupportsFloat | typing.SupportsIndex
    ) -> RefPart:
        """
        Set reference particle kinetic energy [MeV]
        """
    def set_mass_MeV(
        self, mass_MeV: typing.SupportsFloat | typing.SupportsIndex
    ) -> RefPart:
        """
        Set reference particle rest mass * c^2, expressed as an energy [MeV]
        """
    def set_species(self, species_name: str) -> RefPart:
        """
        Set reference particle species by name.

        Sets charge, mass, and gyromagnetic anomaly for a known particle species.
        Returns self for chaining, e.g.: ref.set_species("electron").set_kin_energy_MeV(2.0e3)
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
        reference charge [C]
        """
    @charge.setter
    def charge(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
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
    def gyromagnetic_anomaly(self) -> float:
        """
        reference particle gyromagnetic anomaly [unitless]
        """
    @gyromagnetic_anomaly.setter
    def gyromagnetic_anomaly(
        self, arg0: typing.SupportsFloat | typing.SupportsIndex
    ) -> None: ...
    @property
    def kin_energy_MeV(self) -> float:
        """
        Get reference particle energy [MeV]
        """
    @property
    def mass(self) -> float:
        """
        reference rest mass [kg]
        """
    @mass.setter
    def mass(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def mass_MeV(self) -> float:
        """
        Get reference particle rest mass * c^2, expressed as an energy [MeV]
        """
    @property
    def pt(self) -> float:
        """
        energy deviation, normalized by rest energy
        """
    @pt.setter
    def pt(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def px(self) -> float:
        """
        momentum in x divided by m*c = beta_x*gamma [unitless]
        """
    @px.setter
    def px(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def py(self) -> float:
        """
        momentum in y divided by m*c = beta_y*gamma [unitless]
        """
    @py.setter
    def py(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def pz(self) -> float:
        """
        momentum in z divided by m*c = beta_z*gamma [unitless]
        """
    @pz.setter
    def pz(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def qm_ratio_SI(self) -> float:
        """
        Get reference particle charge to mass ratio [C/kg]
        """
    @property
    def rigidity_Tm(self) -> float:
        """
        Get reference particle magnetic rigidity Brho [T*m]
        """
    @property
    def s(self) -> float:
        """
        integrated orbit path length [m]
        """
    @s.setter
    def s(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def sedge(self) -> float:
        """
        value of s at entrance of the current beamline element
        """
    @sedge.setter
    def sedge(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def t(self) -> float:
        """
        clock time * c [m]
        """
    @t.setter
    def t(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def x(self) -> float:
        """
        horizontal position x [m]
        """
    @x.setter
    def x(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def y(self) -> float:
        """
        vertical position y [m]
        """
    @y.setter
    def y(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def z(self) -> float:
        """
        longitudinal position y [m]
        """
    @z.setter
    def z(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class Envelope:
    envelope: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        arg1: typing.SupportsFloat | typing.SupportsIndex,
    ) -> None: ...
    def beam_moments(self, ref: RefPart) -> dict[str, float]:
        """
        Calculate beam moments (position and momentum moments, emittances,
        Twiss functions, dispersion, ...) from this envelope's covariance
        matrix and a reference particle. The envelope counterpart of
        ImpactXParticleContainer.beam_moments().
        """
    @property
    def beam_intensity(self) -> float: ...
    @beam_intensity.setter
    def beam_intensity(
        self, arg1: typing.SupportsFloat | typing.SupportsIndex
    ) -> Envelope: ...

class CoordSystem:
    """
    Members:

      s

      t
    """

    __members__: typing.ClassVar[dict[str, CoordSystem]]
    s: typing.ClassVar[CoordSystem]
    t: typing.ClassVar[CoordSystem]
    @typing.overload
    def __eq__(self, other: CoordSystem) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: CoordSystem) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class ImpactXParIter(amrex.space3d.amrex_3d_pybind.ParIter_pureSoA_11_0_polymorphic):
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_11_0_polymorphic,
        level: typing.SupportsInt | typing.SupportsIndex,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_11_0_polymorphic,
        level: typing.SupportsInt | typing.SupportsIndex,
        info: amrex.space3d.amrex_3d_pybind.MFItInfo,
    ) -> None: ...
    def pc(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_11_0_polymorphic: ...

class ImpactXParConstIter(
    amrex.space3d.amrex_3d_pybind.ParConstIter_pureSoA_11_0_polymorphic
):
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_11_0_polymorphic,
        level: typing.SupportsInt | typing.SupportsIndex,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        particle_container: amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_11_0_polymorphic,
        level: typing.SupportsInt | typing.SupportsIndex,
        info: amrex.space3d.amrex_3d_pybind.MFItInfo,
    ) -> None: ...
    def pc(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_11_0_polymorphic: ...

class ImpactXParticleContainer(
    amrex.space3d.amrex_3d_pybind.ParticleContainer_pureSoA_11_0_polymorphic
):
    ConstIterator = ImpactXParConstIter
    Iterator = ImpactXParIter
    @staticmethod
    def at_fixed_t(*args, **kwds):
        """
        Temporarily represent the beam at fixed t, e.g., to exchange particles with a code
        that uses time as the independent variable.

        On entry, the particle coordinates are transformed from fixed s to fixed t.
        When the ``with`` block is left, normally or through an exception, they are
        transformed back to fixed s.

        Inside the block, the particle arrays hold ``x, y, z, px, py, pz`` instead of
        ``x, y, t, px, py, pt``: ``z`` is the longitudinal position relative to the
        reference particle, in meters, ``px, py`` are the transverse momenta and ``pz``
        is the deviation from the reference momentum, all normalized by the reference
        momentum.
        The arrays are named accordingly: ``position_z`` and ``momentum_z`` replace
        ``position_t`` and ``momentum_t``.
        The arguments of ``add_n_particles`` keep their names ``t`` and ``pt``; inside the
        block, they take ``z`` and ``pz``.

        Both transformations take the design energy from ``self.ref.pt`` at the moment
        they run: the one on entry uses the reference particle as it is when the block
        starts, the one on exit uses the reference particle as it is when the block ends.
        If the block changes the reference energy, it must also write the particle
        coordinates relative to the new reference particle before the block ends.
        Only ``ref.pt`` is read, so update it (e.g., with ``ref.set_kin_energy_MeV``)
        rather than ``ref.pz`` alone.

        Parameters
        ----------
        self : ImpactXParticleContainer
            The particle container; its particles must be at fixed s.

        Yields
        ------
        ImpactXParticleContainer
            The same particle container, now at fixed t.

        Raises
        ------
        RuntimeError
            If the particles are not at fixed s on entry, or no longer at fixed t
            when the block ends without an exception.

        Examples
        --------
        >>> with sim.beam.at_fixed_t() as beam:
        ...     beam.add_n_particles(x, y, z, px, py, pz, qm, bunch_charge=charge_C)
        """
    def _add_n_particles(
        self,
        x: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        y: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        t: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        px: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        py: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        pt: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
        qm: typing.SupportsFloat | typing.SupportsIndex,
        bunch_charge: typing.SupportsFloat | typing.SupportsIndex | None = None,
        w: amrex.space3d.amrex_3d_pybind.PODVector_real_std | None = None,
        sx: amrex.space3d.amrex_3d_pybind.PODVector_real_std | None = None,
        sy: amrex.space3d.amrex_3d_pybind.PODVector_real_std | None = None,
        sz: amrex.space3d.amrex_3d_pybind.PODVector_real_std | None = None,
    ) -> None:
        """
        Add new particles to the container for fixed s.

        Either the total charge (bunch_charge) or the weight of each
        particle (w) must be provided.
        With bunch_charge, each particle represents the same number of real
        particles with the charge of the reference particle, which must be
        set up before.

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
        :param bunch_charge: magnitude of the total charge of the particles added in this call in C
        :param w: weight of each particle: how many real particles to represent
        :param sx: spin component in x
        :param sy: spin component in y
        :param sz: spin component in z
        """
    def add_n_particles(
        self,
        x,
        y,
        t,
        px,
        py,
        pt,
        qm,
        bunch_charge=None,
        w=None,
        sx=None,
        sy=None,
        sz=None,
    ):
        """
        Add new particles to the container.

        The coordinate and weight arguments accept NumPy or CuPy arrays (or
        array-likes), as well as pyAMReX ``PODVector`` objects. Inputs are copied
        into device-compatible PODVectors as needed.

        Either the total charge (``bunch_charge``) or the weight of each particle
        (``w``) must be provided.
        With ``bunch_charge``, each particle represents the same number of real
        particles with the charge of the reference particle, which must be set up
        before.

        Note: This can only be used *after* the grids have been created, i.e. after
        ``ImpactX.init_grids`` has been called.

        Parameters
        ----------
        x, y, t, px, py, pt : array_like
            Particle coordinates, relative to the reference particle.
            ``x`` and ``y`` are positions in m, ``px`` and ``py`` are momenta
            normalized by the reference momentum (unitless).
            At fixed s (the default), ``t`` is the time-of-flight c*t in m and
            ``pt`` is the energy deviation p_t = -Delta(gamma)/(beta_0 gamma_0)
            (unitless).
            Inside a ``beam.at_fixed_t()`` block, ``t`` takes the position z in m
            and ``pt`` takes the momentum deviation
            p_z = Delta(beta_z gamma)/(beta_0 gamma_0) (unitless).
        qm : float
            Charge over mass in 1/eV.
        bunch_charge : float, optional
            Magnitude of the total charge of the particles added in this call in C.
        w : array_like, optional
            Weight of each particle: how many real particles to represent.
        sx, sy, sz : array_like, optional
            Spin components in x, y, z.
        """
    def beam_moments(self) -> dict[str, float]:
        """
        Calculate beam moments at current ``s`` like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters.
        """
    def beam_moments_history(self):
        """
        Return the history of the beam as calculated by the reduced beam characteristics on every step.
        """
    def beam_moments_history_list(self) -> list[dict[str, float]]:
        """
        Return the history of the beam moments on every step.
        """
    def clear(self, keep_mass: bool = False, keep_charge: bool = False) -> None:
        """
        Empty the container and reset the reference particle
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
    def record_beam_moments(self) -> None:
        """
        Calculate & record the beam moments at current s
        """
    def reduced_beam_characteristics(self) -> dict[str, float]:
        """
        Compute reduced beam characteristics like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters.
        """
    def ref_particle(self) -> RefPart:
        """
        Access the reference particle.

        Deprecated: use ``beam.ref``.
        """
    def reset_beam_moments_history(self) -> None:
        """
        Reset the history of the beam moments.
        """
    def set_bucket_length(
        self, bucket_length: typing.SupportsFloat | typing.SupportsIndex
    ) -> None:
        """
        Set bucket length for particle boundary condition.
        """
    def set_ref_particle(self, refpart: RefPart) -> None:
        """
        Set reference particle attributes.
        """
    @property
    def coord_system(self) -> CoordSystem:
        """
        Get the current coordinate system of particles in this container
        """
    @property
    def ref(self) -> RefPart:
        """
        Access the reference particle.
        """
    @property
    def store_beam_moments(self) -> bool:
        """
        In situ calculate and store the beam moments for every simulation step.
        """
    @store_beam_moments.setter
    def store_beam_moments(self, arg1: bool) -> None: ...

class ImpactX:
    def DistributionMap(
        self, lev: typing.SupportsInt | typing.SupportsIndex
    ) -> amrex.space3d.amrex_3d_pybind.DistributionMapping: ...
    def Geom(
        self, lev: typing.SupportsInt | typing.SupportsIndex
    ) -> amrex.space3d.amrex_3d_pybind.Geometry: ...
    def __init__(self) -> None: ...
    def add_particles(
        self,
        bunch_charge: typing.SupportsFloat | typing.SupportsIndex,
        distr: distribution.Empty
        | distribution.Gaussian
        | distribution.Kurth4D
        | distribution.Kurth6D
        | distribution.KVdist
        | distribution.Thermal
        | distribution.Triangle
        | distribution.Semigaussian
        | distribution.Waterbag,
        npart: typing.SupportsInt | typing.SupportsIndex,
        spin_distr: distribution.SpinvMF | None = None,
    ) -> None:
        """
        Particle tracking mode:Generate and add n particles to the particle container.

        Will also resize the geometry based on the updated particle
        distribution's extent and then redistribute particles in according
        AMReX grid boxes.
        """
    def boxArray(
        self, lev: typing.SupportsInt | typing.SupportsIndex
    ) -> amrex.space3d.amrex_3d_pybind.BoxArray: ...
    def deposit_charge(self) -> None:
        """
        Deposit charge in x,y,z.
        """
    def evolve(self) -> None:
        """
        Run the main simulation loop.
        """
    def finalize(self) -> None:
        """
        Deallocate all contexts and data.
        """
    def init_beam_distribution_from_inputs(self) -> None: ...
    def init_envelope(
        self,
        ref: RefPart,
        distr: distribution.Empty
        | distribution.Gaussian
        | distribution.Kurth4D
        | distribution.Kurth6D
        | distribution.KVdist
        | distribution.Thermal
        | distribution.Triangle
        | distribution.Semigaussian
        | distribution.Waterbag,
        intensity: typing.SupportsFloat | typing.SupportsIndex | None = None,
    ) -> None:
        """
        Envelope tracking mode:Create a 6x6 covariance matrix from a distribution and then initialize the simulation for envelope tracking relative to a reference particle.
        """
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

        Deprecated: use ``sim.beam``.
        """
    def phi(
        self, lev: typing.SupportsInt | typing.SupportsIndex
    ) -> amrex.space3d.amrex_3d_pybind.MultiFab:
        """
        scalar potential per level
        """
    def resize_mesh(self) -> None:
        """
        Resize the mesh :py:attr:`~domain` based on the :py:attr:`~dynamic_size` and related parameters.
        """
    def rho(self, lev: typing.SupportsInt | typing.SupportsIndex) -> typing.Any:
        """
        charge density per level.

        For the 2D space-charge models (``2D`` / ``2.5D``) this returns the
        x-y projected (2D) charge density that is actually solved; for the
        3D models it returns the 3D charge density.
        """
    def space_charge_field(
        self, lev: typing.SupportsInt | typing.SupportsIndex, comp: str
    ) -> amrex.space3d.amrex_3d_pybind.MultiFab:
        """
        space charge force (vector: x,y,z) per level
        """
    def track_envelope(self) -> None:
        """
        Run the envelope tracking simulation loop.
        """
    def track_particles(self) -> None:
        """
        Run the particle tracking simulation loop.
        """
    def track_reference(self, arg0: RefPart) -> None:
        """
        Run the reference orbit tracking simulation loop.
        """
    @property
    def abort_on_unused_inputs(self) -> int:
        """
        Configure simulation to abort AFTER it has run
        if there are unused parameters in the input.
        """
    @abort_on_unused_inputs.setter
    def abort_on_unused_inputs(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
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
    def always_warn_immediately(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def beam(self) -> ImpactXParticleContainer:
        """
        Access the beam particle container.
        """
    @property
    def blocking_factor_x(self) -> list[int]:
        """
        AMReX blocking factor for a direction, per MR level. Also sets max_grid_size_x/_y/_z in this direction to the same value; assign that property afterwards to overwrite. The all-direction max_grid_size does not overwrite it, since per-direction values always take precedence.
        """
    @blocking_factor_x.setter
    def blocking_factor_x(
        self, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]
    ) -> None: ...
    @property
    def blocking_factor_y(self) -> list[int]:
        """
        AMReX blocking factor for a direction, per MR level. Also sets max_grid_size_x/_y/_z in this direction to the same value; assign that property afterwards to overwrite. The all-direction max_grid_size does not overwrite it, since per-direction values always take precedence.
        """
    @blocking_factor_y.setter
    def blocking_factor_y(
        self, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]
    ) -> None: ...
    @property
    def blocking_factor_z(self) -> list[int]:
        """
        AMReX blocking factor for a direction, per MR level. Also sets max_grid_size_x/_y/_z in this direction to the same value; assign that property afterwards to overwrite. The all-direction max_grid_size does not overwrite it, since per-direction values always take precedence.
        """
    @blocking_factor_z.setter
    def blocking_factor_z(
        self, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]
    ) -> None: ...
    @property
    def csr(self) -> bool:
        """
        Enable or disable Coherent Synchrotron Radiation (CSR) calculations (default: disabled).
        """
    @csr.setter
    def csr(self, arg1: bool) -> None: ...
    @property
    def csr_bins(self) -> int:
        """
        Number of longitudinal bins used for CSR calculations (default: 150).
        """
    @csr_bins.setter
    def csr_bins(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def diag_file_min_digits(self) -> int:
        """
        The minimum number of digits (default: 6) used for the step
        number appended to the diagnostic file names.
        """
    @diag_file_min_digits.setter
    def diag_file_min_digits(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def diag_file_prefix(self) -> str:
        """
        Root directory for diagnostic output (default: ``diags``).
        """
    @diag_file_prefix.setter
    def diag_file_prefix(self, arg1: str) -> None: ...
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
    def eigenemittances(self) -> bool:
        """
        Enable or disable eigenemittance diagnostic calculations (default: disabled).
        """
    @eigenemittances.setter
    def eigenemittances(self, arg1: bool) -> None: ...
    @property
    def envelope(self) -> Envelope:
        """
        Access the beam envelope (6x6 covariance matrix and beam intensity)
        used for envelope tracking. Only available after init_envelope().
        """
    @property
    def finest_level(self) -> int:
        """
        The currently finest level of mesh-refinement used. This is always less or equal to max_level.
        """
    @property
    def hook(self) -> UnorderedMap:
        """
        User-defined function hooks that are called, e.g, during tracking.
        """
    @hook.setter
    def hook(self, arg0: UnorderedMap) -> None: ...
    @property
    def isr(self) -> bool:
        """
        Enable or disable Incoherent Synchrotron Radiation (ISR) calculations (default: disabled).
        """
    @isr.setter
    def isr(self, arg1: bool) -> None: ...
    @property
    def isr_on_ref_part(self) -> bool:
        """
        Flag to determine whether ISR radiation loss is applied to the reference particle (default: False).
        """
    @isr_on_ref_part.setter
    def isr_on_ref_part(self, arg1: bool) -> None: ...
    @property
    def isr_order(self) -> int:
        """
        Number of terms in the Taylor series retained for quantum effects (default: 1).
        """
    @isr_order.setter
    def isr_order(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def lattice(self) -> typing.Any:
        """
        Access the accelerator element lattice.
        """
    @lattice.setter
    def lattice(self, arg1: collections.abc.Iterable) -> None: ...
    @property
    def max_grid_size(self) -> list[int]:
        """
        AMReX maximum box size along all directions, per MR level. Boxes larger than this are chopped; aim for one box per parallel process. Per-direction values (max_grid_size_x/_y/_z, also set by blocking_factor_x/_y/_z) take precedence, no matter in which order they are assigned; assign those to change a direction that is already set.
        """
    @max_grid_size.setter
    def max_grid_size(
        self, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]
    ) -> None: ...
    @property
    def max_grid_size_x(self) -> list[int]:
        """
        AMReX maximum box size for a direction, per MR level. Boxes larger than this are chopped, aim for one box per parallel process. Takes precedence over the all-direction max_grid_size, no matter in which order they are assigned.
        """
    @max_grid_size_x.setter
    def max_grid_size_x(
        self, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]
    ) -> None: ...
    @property
    def max_grid_size_y(self) -> list[int]:
        """
        AMReX maximum box size for a direction, per MR level. Boxes larger than this are chopped, aim for one box per parallel process. Takes precedence over the all-direction max_grid_size, no matter in which order they are assigned.
        """
    @max_grid_size_y.setter
    def max_grid_size_y(
        self, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]
    ) -> None: ...
    @property
    def max_grid_size_z(self) -> list[int]:
        """
        AMReX maximum box size for a direction, per MR level. Boxes larger than this are chopped, aim for one box per parallel process. Takes precedence over the all-direction max_grid_size, no matter in which order they are assigned.
        """
    @max_grid_size_z.setter
    def max_grid_size_z(
        self, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]
    ) -> None: ...
    @property
    def max_level(self) -> int:
        """
        The maximum mesh-refinement level for the simulation.
        """
    @max_level.setter
    def max_level(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def memory_profiler(self) -> bool:
        """
        This parameter can be used to disable tiny profiler's memory arena profiling at runtime. If tiny_profiler.enabled is false, this parameter has no effects.
        """
    @memory_profiler.setter
    def memory_profiler(self, arg1: bool) -> None: ...
    @property
    def mlmg_absolute_tolerance(self) -> float:
        """
        The absolute tolerance with which the space-charge fields should be calculated in units of V/m^2. More specifically, the acceptable residual with which the solution can be considered converged. In general this should be left as the default, but in cases where the simulation state changes very little between steps it can occur that the initial guess for the MLMG solver is so close to the converged value that it fails to improve that solution sufficiently to reach the mlmg_relative_tolerance value.
        """
    @mlmg_absolute_tolerance.setter
    def mlmg_absolute_tolerance(
        self, arg1: typing.SupportsFloat | typing.SupportsIndex
    ) -> None: ...
    @property
    def mlmg_max_iters(self) -> int:
        """
        Maximum number of iterations used for MLMG solver for space-charge fields calculation. In case if MLMG converges but fails to reach the desired self_fields_required_precision, this parameter may be increased.
        """
    @mlmg_max_iters.setter
    def mlmg_max_iters(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def mlmg_relative_tolerance(self) -> float:
        """
        The relative precision with which the electrostatic space-charge fields should be calculated. More specifically, the space-charge fields are computed with an iterative Multi-Level Multi-Grid (MLMG) solver. This solver can fail to reach the default precision within a reasonable time.
        """
    @mlmg_relative_tolerance.setter
    def mlmg_relative_tolerance(
        self, arg1: typing.SupportsFloat | typing.SupportsIndex
    ) -> None: ...
    @property
    def mlmg_verbosity(self) -> int:
        """
        The verbosity used for MLMG solver for space-charge fields calculation. Currently MLMG solver looks for verbosity levels from 0-5. A higher number results in more verbose output.
        """
    @mlmg_verbosity.setter
    def mlmg_verbosity(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def n_cell(self) -> list[int]:
        """
        The number of grid points along each direction on the coarsest level.
        """
    @n_cell.setter
    def n_cell(
        self,
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex],
            "FixedSize(3)",
        ],
    ) -> None: ...
    @property
    def omp_threads(self) -> str:
        """
        Controls the number of OpenMP threads to use (ImpactX default: "nosmt").
        https://amrex-codes.github.io/amrex/docs_html/InputsComputeBackends.html.
        """
    @omp_threads.setter
    def omp_threads(
        self, arg1: typing.SupportsInt | typing.SupportsIndex | str
    ) -> None: ...
    @property
    def particle_bc(self) -> str:
        """
        Optional methods to apply a longitudinal particle boundary condition.
        """
    @particle_bc.setter
    def particle_bc(self, arg1: str) -> None: ...
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
    def particle_shape(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def periods(self) -> int:
        """
        The number of periods to repeat the lattice.
        """
    @periods.setter
    def periods(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def poisson_solver(self) -> str:
        """
        The numerical solver to solve the Poisson equation when calculating space charge effects. Either multigrid (default) or fft.
        """
    @poisson_solver.setter
    def poisson_solver(self, arg1: str) -> None: ...
    @property
    def prob_relative(self) -> list[float]:
        """
        The field mesh spans, per direction, multiple times the maximum physical extent of beam particles, as given by this factor.
        """
    @prob_relative.setter
    def prob_relative(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
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
    def space_charge(self) -> str:
        """
        The model to be used when calculating space charge effects. Either off, 2D, 3D, Gauss3D, Gauss2p5D, or 2p5D.
        """
    @space_charge.setter
    def space_charge(self, arg1: bool | str) -> None: ...
    @property
    def space_charge_apply_longitudinal_kick(self) -> bool:
        """
        Enable or disable longitudinal space charge kick in 2.5D space charge solver (default: enabled).
        """
    @space_charge_apply_longitudinal_kick.setter
    def space_charge_apply_longitudinal_kick(self, arg1: bool) -> None: ...
    @property
    def space_charge_gauss_charge_z_bins(self) -> int:
        """
        Number of longitudinal bins (grid points) for computing the linear charge density, at least 2 (default: ``129``).
        """
    @space_charge_gauss_charge_z_bins.setter
    def space_charge_gauss_charge_z_bins(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def space_charge_gauss_long_scale(self) -> float:
        """
        Longitudinal space charge scale for the Gauss2p5D space charge model. Approximation affecting only the longitudinal momentum (``pt``) kick. If not set, it defaults to ``1.103 * gamma * sigma_z``, estimated in-situ from the current reduced beam characteristics, which is a typical value when comparing to a 3D model.
        """
    @space_charge_gauss_long_scale.setter
    def space_charge_gauss_long_scale(
        self, arg1: typing.SupportsFloat | typing.SupportsIndex
    ) -> None: ...
    @property
    def space_charge_gauss_nint(self) -> int:
        """
        Number of steps for computing the integrals (default: ``101``).
        """
    @space_charge_gauss_nint.setter
    def space_charge_gauss_nint(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def space_charge_gauss_taylor_delta(self) -> float:
        """
        Initial region for computing the integrals (default: ``0.01``).
        """
    @space_charge_gauss_taylor_delta.setter
    def space_charge_gauss_taylor_delta(
        self, arg1: typing.SupportsFloat | typing.SupportsIndex
    ) -> None: ...
    @property
    def space_charge_num_longitudinal_bins(self) -> int:
        """
        Number of longitudinal bins (grid points) for 2.5D space charge calculation, at least 2 (default: ``100``).
        """
    @space_charge_num_longitudinal_bins.setter
    def space_charge_num_longitudinal_bins(
        self, arg1: typing.SupportsInt | typing.SupportsIndex
    ) -> None: ...
    @property
    def spin(self) -> bool:
        """
        Enable or disable particle spin tracking (default: disabled).
        """
    @spin.setter
    def spin(self, arg1: bool) -> None: ...
    @property
    def strang_split(self) -> bool:
        """
        Compose the collective effect kicks and the element transport in a second-order, time-symmetric Strang split (default: enabled) or to first order.
        """
    @strang_split.setter
    def strang_split(self, arg1: bool) -> None: ...
    @property
    def tiny_profiler(self) -> bool:
        """
        This parameter can be used to disable tiny profiling including CArena memory profiling at runtime.
        """
    @tiny_profiler.setter
    def tiny_profiler(self, arg1: bool) -> None: ...
    @property
    def tiny_profiler_file(self) -> str:
        """
        If this parameter is empty, the output of tiny profiling is dumped on the default out stream of AMReX. If it's not empty, it specifies the file name for the output. Note that /dev/null is a special name that mean a null file.
        """
    @tiny_profiler_file.setter
    def tiny_profiler_file(self, arg1: str) -> None: ...
    @property
    def tracking_element(
        self,
    ) -> (
        elements.Empty
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
        | elements.ExactCFbend
        | elements.ExactDrift
        | elements.ExactMultipole
        | elements.ExactQuad
        | elements.ExactSbend
        | elements.Kicker
        | elements.LinearMap
        | elements.Marker
        | elements.Multipole
        | elements.NonlinearLens
        | elements.PlaneXYRot
        | elements.PolygonAperture
        | elements.Programmable
        | elements.PRot
        | elements.Quad
        | elements.QuadEdge
        | elements.RFCavity
        | elements.Sbend
        | elements.ShortRF
        | elements.SoftSolenoid
        | elements.SoftQuadrupole
        | elements.Sol
        | elements.Source
        | elements.SpinMap
        | elements.TaperedPL
        | elements.ThinDipole
        | None
    ):
        """
        For tracking hooks/callbacks, the current lattice element.
        """
    @property
    def tracking_period(self) -> int:
        """
        For tracking hooks/callbacks, the period in the lattice (e.g., turn or channel period)
        """
    @property
    def tracking_step(self) -> int:
        """
        For tracking hooks/callbacks, a global step of the simulation.

        A state of internal simulation steps, increments also for space charge slice steps in elements.
        We start in 'step 0' (initial state).
        """
    @property
    def verbose(self) -> int:
        """
        Controls how much information is printed to the terminal, when running ImpactX.
        ``0`` for silent, higher is more verbose. Default is ``1``.
        """
    @verbose.setter
    def verbose(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class UnorderedMap:
    def __bool__(self) -> bool:
        """
        Check whether the map is nonempty
        """
    @typing.overload
    def __contains__(self, arg0: str) -> bool: ...
    @typing.overload
    def __contains__(self, arg0: typing.Any) -> bool: ...
    def __delitem__(self, arg0: str) -> None: ...
    def __getitem__(self, arg0: str) -> collections.abc.Callable[[ImpactX], None]: ...
    def __init__(self) -> None: ...
    def __iter__(self) -> collections.abc.Iterator[str]: ...
    def __len__(self) -> int: ...
    def __setitem__(
        self, arg0: str, arg1: collections.abc.Callable[[ImpactX], None]
    ) -> None: ...
    def items(self) -> typing.ItemsView: ...
    def keys(self) -> typing.KeysView: ...
    def values(self) -> typing.ValuesView: ...

class Config:
    gpu_backend = None
    """
    GPU backend ('CUDA', 'HIP' or 'SYCL'), None without GPU support
    """
    have_fft: typing.ClassVar[bool] = True
    """
    Build supports FFT-based solvers
    """
    have_gpu: typing.ClassVar[bool] = False
    """
    Build supports GPUs
    """
    have_mpi: typing.ClassVar[bool] = True
    """
    Build supports MPI
    """
    have_omp: typing.ClassVar[bool] = True
    """
    Build supports OpenMP
    """
    have_openpmd: typing.ClassVar[bool] = True
    """
    Build supports openPMD I/O
    """
    have_simd: typing.ClassVar[bool] = False
    """
    Build supports explicit SIMD vectorization
    """
    openpmd_backends: typing.ClassVar[dict] = {
        "adios1": False,
        "adios2": False,
        "hdf5": True,
        "json": True,
        "mpi": True,
        "toml": True,
    }
    """
    openPMD I/O backends available in this build
    """
    precision: typing.ClassVar[str] = "DOUBLE"
    """
    Floating point precision of amrex::Real ('SINGLE' or 'DOUBLE')
    """
    precision_particles: typing.ClassVar[str] = "DOUBLE"
    """
    Floating point precision of amrex::ParticleReal ('SINGLE' or 'DOUBLE')
    """
    simd_size: typing.ClassVar[int] = 1
    """
    Number of amrex::ParticleReal elements in a native SIMD vector
    """
    @staticmethod
    def to_dict() -> dict:
        """
        Return the ImpactX build configuration as a dictionary.
        """

def coordinate_transformation(
    pc: ImpactXParticleContainer, direction: CoordSystem
) -> None:
    """
    Transform coordinates from fixed s to fixed t or vice versa, using the reference energy pc.ref.pt at the time of the call.
    """

def create_envelope(
    arg0: distribution.Empty
    | distribution.Gaussian
    | distribution.Kurth4D
    | distribution.Kurth6D
    | distribution.KVdist
    | distribution.Thermal
    | distribution.Triangle
    | distribution.Semigaussian
    | distribution.Waterbag,
    arg1: typing.SupportsFloat | typing.SupportsIndex | None,
) -> Envelope: ...
def flatten_charge_to_2D(
    arg0: ImpactX,
) -> dict[
    int,
    tuple[
        amrex.space3d.amrex_3d_pybind.MultiFab, amrex.space3d.amrex_3d_pybind.MultiFab
    ],
]: ...
@typing.overload
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
    | elements.ExactCFbend
    | elements.ExactDrift
    | elements.ExactMultipole
    | elements.ExactQuad
    | elements.ExactSbend
    | elements.Kicker
    | elements.LinearMap
    | elements.Marker
    | elements.Multipole
    | elements.NonlinearLens
    | elements.PlaneXYRot
    | elements.PolygonAperture
    | elements.Programmable
    | elements.PRot
    | elements.Quad
    | elements.QuadEdge
    | elements.RFCavity
    | elements.Sbend
    | elements.ShortRF
    | elements.SoftSolenoid
    | elements.SoftQuadrupole
    | elements.Sol
    | elements.Source
    | elements.SpinMap
    | elements.TaperedPL
    | elements.ThinDipole,
    step: typing.SupportsInt | typing.SupportsIndex = 0,
    period: typing.SupportsInt | typing.SupportsIndex = 0,
) -> None:
    """
    Push a whole particle beam (incl. reference particle) through an element
    """

@typing.overload
def push(
    ref: RefPart,
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
    | elements.ExactCFbend
    | elements.ExactDrift
    | elements.ExactMultipole
    | elements.ExactQuad
    | elements.ExactSbend
    | elements.Kicker
    | elements.LinearMap
    | elements.Marker
    | elements.Multipole
    | elements.NonlinearLens
    | elements.PlaneXYRot
    | elements.PolygonAperture
    | elements.Programmable
    | elements.PRot
    | elements.Quad
    | elements.QuadEdge
    | elements.RFCavity
    | elements.Sbend
    | elements.ShortRF
    | elements.SoftSolenoid
    | elements.SoftQuadrupole
    | elements.Sol
    | elements.Source
    | elements.SpinMap
    | elements.TaperedPL
    | elements.ThinDipole,
) -> None:
    """
    Push the reference particle through an element
    """

def reverse(
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
    | elements.ExactCFbend
    | elements.ExactDrift
    | elements.ExactMultipole
    | elements.ExactQuad
    | elements.ExactSbend
    | elements.Kicker
    | elements.LinearMap
    | elements.Marker
    | elements.Multipole
    | elements.NonlinearLens
    | elements.PlaneXYRot
    | elements.PolygonAperture
    | elements.Programmable
    | elements.PRot
    | elements.Quad
    | elements.QuadEdge
    | elements.RFCavity
    | elements.Sbend
    | elements.ShortRF
    | elements.SoftSolenoid
    | elements.SoftQuadrupole
    | elements.Sol
    | elements.Source
    | elements.SpinMap
    | elements.TaperedPL
    | elements.ThinDipole,
) -> None:
    """
    Reverse an element in-place so that pushing particles through
    it reverses the effect of the original element.
    """

__author__: str = (
    "Axel Huebl, Chad Mitchell, Ryan Sandberg, Marco Garten, Ji Qiang, et al."
)
__license__: str = "BSD-3-Clause-LBNL"
__version__: str = "26.10"
s: CoordSystem
t: CoordSystem
