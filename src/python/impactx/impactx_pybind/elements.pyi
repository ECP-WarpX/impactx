"""
Accelerator lattice elements in ImpactX
"""

from __future__ import annotations

import typing

import impactx.impactx_pybind

__all__ = [
    "Alignment",
    "Aperture",
    "BeamMonitor",
    "Buncher",
    "CFbend",
    "ChrAcc",
    "ChrDrift",
    "ChrPlasmaLens",
    "ChrQuad",
    "ConstF",
    "DipEdge",
    "Drift",
    "Empty",
    "ExactDrift",
    "ExactSbend",
    "Kicker",
    "KnownElementsList",
    "Multipole",
    "NonlinearLens",
    "PRot",
    "Programmable",
    "Quad",
    "RFCavity",
    "Sbend",
    "ShortRF",
    "SoftQuadrupole",
    "SoftSolenoid",
    "Sol",
    "TaperedPL",
    "Thick",
    "Thin",
    "ThinDipole",
]

class Alignment:
    def __init__(self) -> None:
        """
        Mixin class for lattice elements with horizontal/vertical alignment errors.
        """

    @property
    def dx(self) -> float:
        """
        horizontal translation error in m
        """

    @dx.setter
    def dx(self, arg1: float) -> None: ...
    @property
    def dy(self) -> float:
        """
        vertical translation error in m
        """

    @dy.setter
    def dy(self, arg1: float) -> None: ...
    @property
    def rotation(self) -> float:
        """
        rotation error in the transverse plane in degree
        """

    @rotation.setter
    def rotation(self, arg1: float) -> None: ...

class Aperture(Thin, Alignment):
    def __init__(
        self,
        xmax: float,
        ymax: float,
        shape: str = "rectangular",
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        A short collimator element applying a transverse aperture boundary.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

    @property
    def shape(self) -> str:
        """
        aperture type (rectangular, elliptical)
        """

    @shape.setter
    def shape(self, arg1: str) -> None: ...
    @property
    def xmax(self) -> float:
        """
        maximum horizontal coordinate
        """

    @xmax.setter
    def xmax(self, arg1: float) -> None: ...
    @property
    def ymax(self) -> float:
        """
        maximum vertical coordinate
        """

    @ymax.setter
    def ymax(self, arg1: float) -> None: ...

class BeamMonitor(Thin):
    def __init__(
        self, name: str, backend: str = "default", encoding: str = "g"
    ) -> None:
        """
        This element writes the particle beam out to openPMD data.
        """

    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

    @property
    def alpha(self) -> float:
        """
        Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.
        """

    @alpha.setter
    def alpha(self, arg1: float) -> None: ...
    @property
    def beta(self) -> float:
        """
        Twiss beta (in meters) of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.
        """

    @beta.setter
    def beta(self, arg1: float) -> None: ...
    @property
    def cn(self) -> float:
        """
        Scale factor (in meters^(1/2)) of the IOTA nonlinear magnetic insert element used for computing H and I.
        """

    @cn.setter
    def cn(self, arg1: float) -> None: ...
    @property
    def name(self) -> str:
        """
        name of the series
        """

    @property
    def nonlinear_lens_invariants(self) -> bool:
        """
        Compute and output the invariants H and I within the nonlinear magnetic insert element
        """

    @nonlinear_lens_invariants.setter
    def nonlinear_lens_invariants(self, arg1: bool) -> None: ...
    @property
    def tn(self) -> float:
        """
        Dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I.
        """

    @tn.setter
    def tn(self, arg1: float) -> None: ...

class Buncher(Thin, Alignment):
    def __init__(
        self, V: float, k: float, dx: float = 0, dy: float = 0, rotation: float = 0
    ) -> None:
        """
        A short linear RF cavity element at zero-crossing for bunching.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class CFbend(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        rc: float,
        k: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        An ideal combined function bend (sector bend with quadrupole component).
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class ChrAcc(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        ez: float,
        bz: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        A region of Uniform Acceleration, with chromatic effects included.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

    @property
    def bz(self) -> float:
        """
        magnetic field strength in 1/m
        """

    @bz.setter
    def bz(self, arg1: float) -> None: ...
    @property
    def ez(self) -> float:
        """
        electric field strength in 1/m
        """

    @ez.setter
    def ez(self, arg1: float) -> None: ...

class ChrDrift(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        A Drift with chromatic effects included.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class ChrPlasmaLens(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        k: float,
        units: int = 0,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        An active Plasma Lens with chromatic effects included.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

    @property
    def k(self) -> float:
        """
        focusing strength in 1/m^2 (or T/m)
        """

    @k.setter
    def k(self, arg1: float) -> None: ...
    @property
    def units(self) -> int:
        """
        unit specification for focusing strength
        """

    @units.setter
    def units(self, arg1: int) -> None: ...

class ChrQuad(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        k: float,
        units: int = 0,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        A Quadrupole magnet with chromatic effects included.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

    @property
    def k(self) -> float:
        """
        quadrupole strength in 1/m^2 (or T/m)
        """

    @k.setter
    def k(self, arg1: float) -> None: ...
    @property
    def units(self) -> int:
        """
        unit specification for quad strength
        """

    @units.setter
    def units(self, arg1: int) -> None: ...

class ConstF(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        kx: float,
        ky: float,
        kt: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        A linear Constant Focusing element.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

    @property
    def kt(self) -> float:
        """
        focusing t strength in 1/m
        """

    @kt.setter
    def kt(self, arg1: float) -> None: ...
    @property
    def kx(self) -> float:
        """
        focusing x strength in 1/m
        """

    @kx.setter
    def kx(self, arg1: float) -> None: ...
    @property
    def ky(self) -> float:
        """
        focusing y strength in 1/m
        """

    @ky.setter
    def ky(self, arg1: float) -> None: ...

class DipEdge(Thin, Alignment):
    def __init__(
        self,
        psi: float,
        rc: float,
        g: float,
        K2: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        Edge focusing associated with bend entry or exit.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class Drift(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        A drift.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class Empty(Thin):
    def __init__(self) -> None:
        """
        This element does nothing.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class ExactDrift(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        A Drift using the exact nonlinear map.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class ExactSbend(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        phi: float,
        B: float = 0.0,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        An ideal sector bend using the exact nonlinear map.  When B = 0, the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class Kicker(Thin, Alignment):
    def __init__(
        self,
        xkick: float,
        ykick: float,
        units: str = "dimensionless",
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        A thin transverse kicker element. Kicks are for units "dimensionless" or in "T-m".
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class KnownElementsList:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: (
            Empty
            | Aperture
            | Buncher
            | CFbend
            | ChrAcc
            | ChrDrift
            | ChrPlasmaLens
            | ChrQuad
            | ConstF
            | BeamMonitor
            | DipEdge
            | Drift
            | ExactDrift
            | ExactSbend
            | Kicker
            | Multipole
            | NonlinearLens
            | Programmable
            | PRot
            | Quad
            | RFCavity
            | Sbend
            | ShortRF
            | SoftSolenoid
            | SoftQuadrupole
            | Sol
            | TaperedPL
            | ThinDipole
        ),
    ) -> None: ...
    @typing.overload
    def __init__(self, arg0: list) -> None: ...
    def __iter__(
        self,
    ) -> typing.Iterator[
        Empty
        | Aperture
        | Buncher
        | CFbend
        | ChrAcc
        | ChrDrift
        | ChrPlasmaLens
        | ChrQuad
        | ConstF
        | BeamMonitor
        | DipEdge
        | Drift
        | ExactDrift
        | ExactSbend
        | Kicker
        | Multipole
        | NonlinearLens
        | Programmable
        | PRot
        | Quad
        | RFCavity
        | Sbend
        | ShortRF
        | SoftSolenoid
        | SoftQuadrupole
        | Sol
        | TaperedPL
        | ThinDipole
    ]: ...
    def __len__(self) -> int:
        """
        The length of the list.
        """

    def append(
        self,
        arg0: (
            Empty
            | Aperture
            | Buncher
            | CFbend
            | ChrAcc
            | ChrDrift
            | ChrPlasmaLens
            | ChrQuad
            | ConstF
            | BeamMonitor
            | DipEdge
            | Drift
            | ExactDrift
            | ExactSbend
            | Kicker
            | Multipole
            | NonlinearLens
            | Programmable
            | PRot
            | Quad
            | RFCavity
            | Sbend
            | ShortRF
            | SoftSolenoid
            | SoftQuadrupole
            | Sol
            | TaperedPL
            | ThinDipole
        ),
    ) -> None:
        """
        Add a single element to the list.
        """

    def clear(self) -> None:
        """
        Clear the list to become empty.
        """

    @typing.overload
    def extend(self, arg0: KnownElementsList) -> KnownElementsList:
        """
        Add a list of elements to the list.
        """

    @typing.overload
    def extend(self, arg0: list) -> KnownElementsList:
        """
        Add a list of elements to the list.
        """

    def load_file(self, madx_file, nslice=1): ...
    def pop_back(self) -> None:
        """
        Return and remove the last element of the list.
        """

class Multipole(Thin, Alignment):
    def __init__(
        self,
        multipole: int,
        K_normal: float,
        K_skew: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        A general thin multipole element.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class NonlinearLens(Thin, Alignment):
    def __init__(
        self,
        knll: float,
        cnll: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        Single short segment of the nonlinear magnetic insert element.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class PRot(Thin):
    def __init__(self, phi_in: float, phi_out: float) -> None:
        """
        An exact pole-face rotation in the x-z plane. Both angles are in degrees.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class Programmable:
    ds: float
    nslice: int
    def __init__(self, ds: float = 0.0, nslice: int = 1) -> None:
        """
        A programmable beam optics element.
        """

    def __repr__(self) -> str: ...
    @property
    def beam_particles(
        self,
    ) -> typing.Callable[
        [impactx.impactx_pybind.ImpactXParIter, impactx.impactx_pybind.RefPart], None
    ]:
        """
        hook for beam particles (pti, RefPart)
        """

    @beam_particles.setter
    def beam_particles(
        self,
        arg1: typing.Callable[
            [impactx.impactx_pybind.ImpactXParIter, impactx.impactx_pybind.RefPart],
            None,
        ],
    ) -> None: ...
    @property
    def push(
        self,
    ) -> typing.Callable[[impactx.impactx_pybind.ImpactXParticleContainer, int], None]:
        """
        hook for push of whole container (pc, step)
        """

    @push.setter
    def push(
        self,
        arg1: typing.Callable[
            [impactx.impactx_pybind.ImpactXParticleContainer, int], None
        ],
    ) -> None: ...
    @property
    def ref_particle(self) -> typing.Callable[[impactx.impactx_pybind.RefPart], None]:
        """
        hook for reference particle (RefPart)
        """

    @ref_particle.setter
    def ref_particle(
        self, arg1: typing.Callable[[impactx.impactx_pybind.RefPart], None]
    ) -> None: ...
    @property
    def threadsafe(self) -> bool:
        """
        allow threading via OpenMP for the particle iterator loop, default=False (note: if OMP backend is active)
        """

    @threadsafe.setter
    def threadsafe(self, arg1: bool) -> None: ...

class Quad(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        k: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        A Quadrupole magnet.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class RFCavity(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        escale: float,
        freq: float,
        phase: float,
        cos_coefficients: list[float],
        sin_coefficients: list[float],
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        mapsteps: int = 1,
        nslice: int = 1,
    ) -> None:
        """
        An RF cavity.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class Sbend(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        rc: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        An ideal sector bend.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class ShortRF(Thin, Alignment):
    def __init__(
        self,
        V: float,
        freq: float,
        phase: float = -90.0,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        A short RF cavity element.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class SoftQuadrupole(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        gscale: float,
        cos_coefficients: list[float],
        sin_coefficients: list[float],
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        mapsteps: int = 1,
        nslice: int = 1,
    ) -> None:
        """
        A soft-edge quadrupole.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class SoftSolenoid(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        bscale: float,
        cos_coefficients: list[float],
        sin_coefficients: list[float],
        unit: float = 0,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        mapsteps: int = 1,
        nslice: int = 1,
    ) -> None:
        """
        A soft-edge solenoid.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class Sol(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        ks: float,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
        nslice: int = 1,
    ) -> None:
        """
        An ideal hard-edge Solenoid magnet.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class TaperedPL(Thin, Alignment):
    def __init__(
        self,
        k: float,
        taper: float,
        units: int = 0,
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        A thin nonlinear plasma lens with transverse (horizontal) taper

                     .. math::

                        B_x = g \\left( y + \\frac{xy}{D_x} \\right), \\quad \\quad B_y = -g \\left(x + \\frac{x^2 + y^2}{2 D_x} \\right)

                     where :math:`g` is the (linear) field gradient in T/m and :math:`D_x` is the targeted horizontal dispersion in m.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """

class Thick:
    def __init__(self, ds: float, nslice: float = 1) -> None:
        """
        Mixin class for lattice elements with finite length.
        """

    @property
    def ds(self) -> float:
        """
        segment length in m
        """

    @ds.setter
    def ds(self, arg1: float) -> None: ...
    @property
    def nslice(self) -> int:
        """
        number of slices used for the application of space charge
        """

    @nslice.setter
    def nslice(self, arg1: int) -> None: ...

class Thin:
    def __init__(self) -> None:
        """
        Mixin class for lattice elements with zero length.
        """

    @property
    def ds(self) -> float:
        """
        segment length in m
        """

    @property
    def nslice(self) -> int:
        """
        number of slices used for the application of space charge
        """

class ThinDipole(Thin, Alignment):
    def __init__(
        self, theta: float, rc: float, dx: float = 0, dy: float = 0, rotation: float = 0
    ) -> None:
        """
        A thin kick model of a dipole bend.
        """

    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
