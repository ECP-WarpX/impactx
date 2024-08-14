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
    @property
    def V(self) -> float:
        """
        Normalized RF voltage drop V = Emax*L/(c*Brho)
        """
    @V.setter
    def V(self, arg1: float) -> None: ...
    @property
    def k(self) -> float:
        """
        Wavenumber of RF in 1/m
        """
    @k.setter
    def k(self, arg1: float) -> None: ...

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
    @property
    def k(self) -> float:
        """
        Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing
        """
    @k.setter
    def k(self, arg1: float) -> None: ...
    @property
    def rc(self) -> float:
        """
        Radius of curvature in m
        """
    @rc.setter
    def rc(self, arg1: float) -> None: ...

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
        unit: int = 0,
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
    def unit(self) -> int:
        """
        unit specification for focusing strength
        """
    @unit.setter
    def unit(self, arg1: int) -> None: ...

class ChrQuad(Thick, Alignment):
    def __init__(
        self,
        ds: float,
        k: float,
        unit: int = 0,
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
    def unit(self) -> int:
        """
        unit specification for quad strength
        """
    @unit.setter
    def unit(self, arg1: int) -> None: ...

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
    @property
    def K2(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K2.setter
    def K2(self, arg1: float) -> None: ...
    @property
    def g(self) -> float:
        """
        Gap parameter in m
        """
    @g.setter
    def g(self, arg1: float) -> None: ...
    @property
    def psi(self) -> float:
        """
        Pole face angle in rad
        """
    @psi.setter
    def psi(self, arg1: float) -> None: ...
    @property
    def rc(self) -> float:
        """
        Radius of curvature in m
        """
    @rc.setter
    def rc(self, arg1: float) -> None: ...

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
    @property
    def B(self) -> float:
        """
        Magnetic field in Tesla; when B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B
        """
    @B.setter
    def B(self, arg1: float) -> None: ...
    @property
    def phi(self) -> float:
        """
        Bend angle in degrees
        """
    @phi.setter
    def phi(self, arg1: float) -> None: ...

class Kicker(Thin, Alignment):
    def __init__(
        self,
        xkick: float,
        ykick: float,
        unit: str = "dimensionless",
        dx: float = 0,
        dy: float = 0,
        rotation: float = 0,
    ) -> None:
        """
        A thin transverse kicker element. Kicks are for unit "dimensionless" or in "T-m".
        """
    def __repr__(self) -> str: ...
    def push(
        self, pc: impactx.impactx_pybind.ImpactXParticleContainer, step: int = 0
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @property
    def xkick(self) -> float:
        """
        horizontal kick strength (dimensionless OR T-m)
        """
    @xkick.setter
    def xkick(self, arg1: float) -> None: ...
    @property
    def ykick(self) -> float:
        """
        vertical kick strength (dimensionless OR T-m)
        """
    @ykick.setter
    def ykick(self, arg1: float) -> None: ...

class KnownElementsList:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: Empty
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
        | ThinDipole,
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
        arg0: Empty
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
        | ThinDipole,
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
    @property
    def K_normal(self) -> float:
        """
        Integrated normal multipole coefficient (1/meter^m)
        """
    @K_normal.setter
    def K_normal(self, arg1: float) -> None: ...
    @property
    def K_skew(self) -> float:
        """
        Integrated skew multipole coefficient (1/meter^m)
        """
    @K_skew.setter
    def K_skew(self, arg1: float) -> None: ...
    @property
    def multipole(self) -> int:
        """
        index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
        """
    @multipole.setter
    def multipole(self, arg1: float) -> None: ...

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
    @property
    def cnll(self) -> float:
        """
        distance of singularities from the origin (m)
        """
    @cnll.setter
    def cnll(self, arg1: float) -> None: ...
    @property
    def knll(self) -> float:
        """
        integrated strength of the nonlinear lens (m)
        """
    @knll.setter
    def knll(self, arg1: float) -> None: ...

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
    @property
    def phi_in(self) -> float:
        """
        angle of the reference particle with respect to the longitudinal (z) axis in the original frame in degrees
        """
    @phi_in.setter
    def phi_in(self, arg1: float) -> None: ...
    @property
    def phi_out(self) -> float:
        """
        angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in degrees
        """
    @phi_out.setter
    def phi_out(self, arg1: float) -> None: ...

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
    @property
    def k(self) -> float:
        """
        Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing
        """
    @k.setter
    def k(self, arg1: float) -> None: ...

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
    @property
    def escale(self) -> float:
        """
        scaling factor for on-axis RF electric field in 1/m = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV)
        """
    @escale.setter
    def escale(self, arg1: float) -> None: ...
    @property
    def freq(self) -> float:
        """
        RF frequency in Hz
        """
    @freq.setter
    def freq(self, arg1: float) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: int) -> None: ...
    @property
    def phase(self) -> float:
        """
        RF driven phase in degrees
        """
    @phase.setter
    def phase(self, arg1: float) -> None: ...

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
    @property
    def rc(self) -> float:
        """
        Radius of curvature in m
        """
    @rc.setter
    def rc(self, arg1: float) -> None: ...

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
    @property
    def V(self) -> float:
        """
        Normalized RF voltage V = maximum energy gain/(m*c^2)
        """
    @V.setter
    def V(self, arg1: float) -> None: ...
    @property
    def freq(self) -> float:
        """
        RF frequency in Hz
        """
    @freq.setter
    def freq(self, arg1: float) -> None: ...
    @property
    def phase(self) -> float:
        """
        RF synchronous phase in degrees (phase = 0 corresponds to maximum energy gain, phase = -90 corresponds go zero energy gain for bunching)
        """
    @phase.setter
    def phase(self, arg1: float) -> None: ...

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
    @property
    def gscale(self) -> float:
        """
        Scaling factor for on-axis field gradient in inverse meters
        """
    @gscale.setter
    def gscale(self, arg1: float) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: int) -> None: ...

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
    @property
    def bscale(self) -> float:
        """
        Scaling factor for on-axis magnetic field Bz in inverse meters (if unit = 0) or magnetic field Bz in T (SI units, if unit = 1)
        """
    @bscale.setter
    def bscale(self, arg1: float) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: int) -> None: ...
    @property
    def unit(self) -> int:
        """
        specification of units for scaling of the on-axis longitudinal magnetic field
        """
    @unit.setter
    def unit(self, arg1: float) -> None: ...

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
    @property
    def ks(self) -> float:
        """
        Solenoid strength in m^(-1) (MADX convention) in (magnetic field Bz in T) / (rigidity in T-m)
        """
    @ks.setter
    def ks(self, arg1: float) -> None: ...

class TaperedPL(Thin, Alignment):
    def __init__(
        self,
        k: float,
        taper: float,
        unit: int = 0,
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
    @property
    def k(self) -> float:
        """
        integrated focusing strength in m^(-1) (if unit = 0) or integrated focusing strength in T (if unit = 1)
        """
    @k.setter
    def k(self, arg1: float) -> None: ...
    @property
    def taper(self) -> float:
        """
        horizontal taper parameter in m^(-1) = 1 / (target horizontal dispersion in m)
        """
    @taper.setter
    def taper(self, arg1: float) -> None: ...
    @property
    def unit(self) -> int:
        """
        specification of units for plasma lens focusing strength
        """
    @unit.setter
    def unit(self, arg1: int) -> None: ...

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
    @property
    def rc(self) -> float:
        """
        Effective curvature radius (meters)
        """
    @rc.setter
    def rc(self, arg1: float) -> None: ...
    @property
    def theta(self) -> float:
        """
        Bend angle (degrees)
        """
    @theta.setter
    def theta(self, arg1: float) -> None: ...
