"""
Accelerator lattice elements in ImpactX
"""

from __future__ import annotations

import collections.abc
import typing

import amrex.space3d.amrex_3d_pybind
import impactx.impactx_pybind
from impactx.extensions.Elements import isclose

from . import mixin, transformation

__all__: list[str] = [
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
    "ExactCFbend",
    "ExactDrift",
    "ExactMultipole",
    "ExactQuad",
    "ExactSbend",
    "FilteredElementsList",
    "Kicker",
    "KnownElementsList",
    "LinearMap",
    "Marker",
    "Multipole",
    "NonlinearLens",
    "PRot",
    "PlaneXYRot",
    "PolygonAperture",
    "Programmable",
    "Quad",
    "QuadEdge",
    "RFCavity",
    "Sbend",
    "ShortRF",
    "SoftQuadrupole",
    "SoftSolenoid",
    "Sol",
    "Source",
    "SpinMap",
    "TaperedPL",
    "ThinDipole",
    "isclose",
    "mixin",
    "transformation",
]

class BeamMonitor(mixin.Thin):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        name: str,
        backend: str = "default",
        encoding: str = "g",
        period_sample_intervals: typing.SupportsInt | typing.SupportsIndex = 1,
        particles: bool = True,
        beam_moments: bool = True,
    ) -> None:
        """
        This element writes the particle beam out to openPMD data.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def alpha(self) -> float:
        """
        Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.
        """
    @alpha.setter
    def alpha(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def backend(self) -> str:
        """
        openPMD file backend (e.g. default, bp4, h5)
        """
    @property
    def beam_moments(self) -> bool:
        """
        Output the beam moments (reduced beam characteristics) as openPMD attributes (default: True)
        """
    @beam_moments.setter
    def beam_moments(self, arg1: bool) -> None: ...
    @property
    def beta(self) -> float:
        """
        Twiss beta (in meters) of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.
        """
    @beta.setter
    def beta(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def cn(self) -> float:
        """
        Scale factor (in meters^(1/2)) of the IOTA nonlinear magnetic insert element used for computing H and I.
        """
    @cn.setter
    def cn(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def encoding(self) -> str:
        """
        openPMD iteration encoding: "v" variable-based, "f" file-based, "g" group-based
        """
    @property
    def has_name(self) -> bool: ...
    @property
    def name(self) -> str:
        """
        name of the series
        """
    @property
    def nonlinear_lens_invariants(self) -> bool:
        """
        Compute and output the invariants H and I within the nonlinear magnetic insert element. Requires particles=True.
        """
    @nonlinear_lens_invariants.setter
    def nonlinear_lens_invariants(self, arg1: bool) -> None: ...
    @property
    def particles(self) -> bool:
        """
        Output the individual particles of the beam (default: True)
        """
    @particles.setter
    def particles(self, arg1: bool) -> None: ...
    @property
    def period_sample_intervals(self) -> int:
        """
        for periodic lattices, only output every Nth period (turn or cycle)
        """
    @property
    def tn(self) -> float:
        """
        Dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I.
        """
    @tn.setter
    def tn(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class Aperture(mixin.Named, mixin.Thin, mixin.Alignment):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex,
        repeat_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        repeat_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        shift_odd_x: bool = False,
        shape: str = "rectangular",
        action: str = "transmit",
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A short collimator element applying a transverse aperture boundary.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def action(self) -> str:
        """
        action type (transmit, absorb)
        """
    @action.setter
    def action(self, arg1: str) -> None: ...
    @property
    def aperture_x(self) -> float:
        """
        maximum horizontal coordinate in m; zero or less removes the horizontal boundary
        """
    @aperture_x.setter
    def aperture_x(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def aperture_y(self) -> float:
        """
        maximum vertical coordinate in m; zero or less removes the vertical boundary
        """
    @aperture_y.setter
    def aperture_y(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def repeat_x(self) -> float:
        """
        horizontal period for repeated aperture masking
        """
    @repeat_x.setter
    def repeat_x(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def repeat_y(self) -> float:
        """
        vertical period for repeated aperture masking
        """
    @repeat_y.setter
    def repeat_y(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def shape(self) -> str:
        """
        aperture type (rectangular, elliptical)
        """
    @shape.setter
    def shape(self, arg1: str) -> None: ...
    @property
    def shift_odd_x(self) -> bool:
        """
        for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
        """
    @shift_odd_x.setter
    def shift_odd_x(self, arg1: bool) -> None: ...

class ChrDrift(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A Drift with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """

class ChrQuad(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        k: typing.SupportsFloat | typing.SupportsIndex,
        unit: typing.SupportsInt | typing.SupportsIndex = 0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A Quadrupole magnet with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def k(self) -> float:
        """
        quadrupole strength in 1/m^2 (or T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for quad strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class ChrPlasmaLens(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        k: typing.SupportsFloat | typing.SupportsIndex,
        unit: typing.SupportsInt | typing.SupportsIndex = 0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        An active Plasma Lens with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def k(self) -> float:
        """
        focusing strength in 1/m^2 (or T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for focusing strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class ChrAcc(mixin.Named, mixin.Thick, mixin.Alignment):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        ez: typing.SupportsFloat | typing.SupportsIndex,
        bz: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A region of Uniform Acceleration, with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def bz(self) -> float:
        """
        magnetic field strength in 1/m
        """
    @bz.setter
    def bz(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def ez(self) -> float:
        """
        electric field strength in 1/m
        """
    @ez.setter
    def ez(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class ConstF(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        kx: typing.SupportsFloat | typing.SupportsIndex,
        ky: typing.SupportsFloat | typing.SupportsIndex,
        kt: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A linear Constant Focusing element.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def kt(self) -> float:
        """
        focusing t strength in 1/m
        """
    @kt.setter
    def kt(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def kx(self) -> float:
        """
        focusing x strength in 1/m
        """
    @kx.setter
    def kx(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def ky(self) -> float:
        """
        focusing y strength in 1/m
        """
    @ky.setter
    def ky(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class DipEdge(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        psi: typing.SupportsFloat | typing.SupportsIndex,
        rc: typing.SupportsFloat | typing.SupportsIndex,
        g: typing.SupportsFloat | typing.SupportsIndex,
        R: typing.SupportsFloat | typing.SupportsIndex = 1.0,
        K0: typing.SupportsFloat | typing.SupportsIndex = 1.6449340668482264,
        K1: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        K2: typing.SupportsFloat | typing.SupportsIndex = 1.0,
        K3: typing.SupportsFloat | typing.SupportsIndex = 0.16666666666666666,
        K4: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        K5: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        K6: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        model: str = "linear",
        location: str = "entry",
        modify_ref_part: bool = False,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        Edge focusing associated with bend entry or exit.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def K0(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K0.setter
    def K0(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def K1(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K1.setter
    def K1(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def K2(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K2.setter
    def K2(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def K3(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K3.setter
    def K3(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def K4(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K4.setter
    def K4(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def K5(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K5.setter
    def K5(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def K6(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K6.setter
    def K6(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def R(self) -> float:
        """
        Length scale for field integrals in m
        """
    @R.setter
    def R(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def g(self) -> float:
        """
        Gap parameter in m
        """
    @g.setter
    def g(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def location(self) -> str:
        """
        Fringe field location (entry or exit)
        """
    @location.setter
    def location(self, arg1: str) -> None: ...
    @property
    def model(self) -> str:
        """
        Fringe field model (linear or nonlinear)
        """
    @model.setter
    def model(self, arg1: str) -> None: ...
    @property
    def modify_ref_part(self) -> bool:
        """
        Apply DipEdge to reference particle (boolean).
        """
    @modify_ref_part.setter
    def modify_ref_part(self, arg1: bool) -> None: ...
    @property
    def psi(self) -> float:
        """
        Pole face angle in rad
        """
    @psi.setter
    def psi(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def rc(self) -> float:
        """
        Radius of curvature in m
        """
    @rc.setter
    def rc(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class QuadEdge(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        k: typing.SupportsFloat | typing.SupportsIndex,
        unit: typing.SupportsInt | typing.SupportsIndex = 0,
        flag: str = "entry",
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A thin quadrupole fringe field element. Flag must be "entry" or "exit".
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def k(self) -> float:
        """
        quadrupole focusing strength (1/meter^2 OR T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for quad strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class Drift(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A drift.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """

class ExactDrift(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A Drift using the exact nonlinear map.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """

class ExactMultipole(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        k_normal: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
        k_skew: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
        unit: typing.SupportsInt | typing.SupportsIndex = 0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        int_order: typing.SupportsInt | typing.SupportsIndex = 2,
        mapsteps: typing.SupportsInt | typing.SupportsIndex = 10,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A thick Multipole magnet using the exact nonlinear Hamiltonian.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def set_coefficients(
        self,
        k_normal: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
        k_skew: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None:
        """
        Replace both arrays at once.

        They are validated together and applied as one step, so a rejected update
        leaves the element unchanged.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def int_order(self) -> int:
        """
        order of symplectic integration used for particle push in applied fields
        """
    @int_order.setter
    def int_order(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def k_normal(self) -> list[float]:
        """
        normal multipole coefficients

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @k_normal.setter
    def k_normal(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def k_skew(self) -> list[float]:
        """
        skew multipole coefficients

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @k_skew.setter
    def k_skew(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for particle push in the applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for multipole strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class ExactCFbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        k_normal: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
        k_skew: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
        unit: typing.SupportsInt | typing.SupportsIndex = 0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        int_order: typing.SupportsInt | typing.SupportsIndex = 2,
        mapsteps: typing.SupportsInt | typing.SupportsIndex = 10,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A thick combined function bending magnet using the exact nonlinear Hamiltonian.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def set_coefficients(
        self,
        k_normal: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
        k_skew: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None:
        """
        Replace both arrays at once.

        They are validated together and applied as one step, so a rejected update
        leaves the element unchanged.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def int_order(self) -> int:
        """
        order of symplectic integration used for particle push in applied fields
        """
    @int_order.setter
    def int_order(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def k_normal(self) -> list[float]:
        """
        normal multipole coefficients

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @k_normal.setter
    def k_normal(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def k_skew(self) -> list[float]:
        """
        skew multipole coefficients

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @k_skew.setter
    def k_skew(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for particle push in the applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for multipole strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class ExactQuad(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        k: typing.SupportsFloat | typing.SupportsIndex,
        unit: typing.SupportsInt | typing.SupportsIndex = 0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        int_order: typing.SupportsInt | typing.SupportsIndex = 2,
        mapsteps: typing.SupportsInt | typing.SupportsIndex = 10,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A Quadrupole magnet using the exact nonlinear Hamiltonian.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def int_order(self) -> int:
        """
        order of symplectic integration used for particle push in applied fields
        """
    @int_order.setter
    def int_order(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def k(self) -> float:
        """
        quadrupole strength in 1/m^2 (or T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for particle push in the applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for quad strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class ExactSbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        phi: typing.SupportsFloat | typing.SupportsIndex,
        B: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal sector bend using the exact nonlinear map.  When B = 0, the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def rc(self, ref: impactx.impactx_pybind.RefPart) -> float:
        """
        Radius of curvature in m
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self, in_degrees: bool = False
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def B(self) -> float:
        """
        Magnetic field in Tesla; when B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B
        """
    @B.setter
    def B(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def phi(self) -> float:
        """
        Bend angle in radian
        """
    @phi.setter
    def phi(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class Kicker(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        xkick: typing.SupportsFloat | typing.SupportsIndex,
        ykick: typing.SupportsFloat | typing.SupportsIndex,
        unit: str = "dimensionless",
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A thin transverse kicker element. Kicks are for unit "dimensionless" or in "T-m".
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def xkick(self) -> float:
        """
        horizontal kick strength (dimensionless OR T-m)
        """
    @xkick.setter
    def xkick(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def ykick(self) -> float:
        """
        vertical kick strength (dimensionless OR T-m)
        """
    @ykick.setter
    def ykick(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class Multipole(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        multipole: typing.SupportsInt | typing.SupportsIndex,
        K_normal: typing.SupportsFloat | typing.SupportsIndex,
        K_skew: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A general thin multipole element.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def K_normal(self) -> float:
        """
        Integrated normal multipole coefficient (1/meter^m)
        """
    @K_normal.setter
    def K_normal(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def K_skew(self) -> float:
        """
        Integrated skew multipole coefficient (1/meter^m)
        """
    @K_skew.setter
    def K_skew(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def multipole(self) -> int:
        """
        index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
        """
    @multipole.setter
    def multipole(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class Empty(mixin.Named, mixin.Thin):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(self) -> None:
        """
        This element does nothing.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """

class Marker(mixin.Named, mixin.Thin):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(self, name: str) -> None:
        """
        This named element does nothing.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """

class NonlinearLens(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        knll: typing.SupportsFloat | typing.SupportsIndex,
        cnll: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        Single short segment of the nonlinear magnetic insert element.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def cnll(self) -> float:
        """
        distance of singularities from the origin (m)
        """
    @cnll.setter
    def cnll(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def knll(self) -> float:
        """
        integrated strength of the nonlinear lens (m)
        """
    @knll.setter
    def knll(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class PlaneXYRot(mixin.Named, mixin.Thin, mixin.Alignment):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        angle: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A rotation in the x-y plane.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self, in_degrees: bool = False
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def angle(self) -> float:
        """
        Rotation angle (rad).
        """
    @angle.setter
    def angle(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class PolygonAperture(mixin.Named, mixin.Thin, mixin.Alignment):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        vertices_x: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
        vertices_y: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
        min_radius2: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        repeat_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        repeat_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        shift_odd_x: bool = False,
        action: str = "transmit",
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A short collimator element described by a polygon with vertices given by their x and y coordinates.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def set_vertices(
        self,
        vertices_x: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
        vertices_y: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
    ) -> None:
        """
        Replace both arrays at once.

        They are validated together and applied as one step, so a rejected update
        leaves the element unchanged.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def action(self) -> str:
        """
        action type (transmit, absorb)
        """
    @action.setter
    def action(self, arg1: str) -> None: ...
    @property
    def min_radius2(self) -> float:
        """
        All particles with radius squared smaller than min_radius2 pass the aperture
        """
    @min_radius2.setter
    def min_radius2(
        self, arg1: typing.SupportsFloat | typing.SupportsIndex
    ) -> None: ...
    @property
    def repeat_x(self) -> float:
        """
        horizontal period for repeated aperture masking
        """
    @repeat_x.setter
    def repeat_x(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def repeat_y(self) -> float:
        """
        vertical period for repeated aperture masking
        """
    @repeat_y.setter
    def repeat_y(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def shift_odd_x(self) -> bool:
        """
        for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
        """
    @shift_odd_x.setter
    def shift_odd_x(self, arg1: bool) -> None: ...
    @property
    def vertices_x(self) -> list[float]:
        """
        horizontal vertex coordinates in m

        Both arrays are always the same length. To change the length, set both at
        once with ``set_vertices()``.
        """
    @vertices_x.setter
    def vertices_x(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def vertices_y(self) -> list[float]:
        """
        vertical vertex coordinates in m

        Both arrays are always the same length. To change the length, set both at
        once with ``set_vertices()``.
        """
    @vertices_y.setter
    def vertices_y(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...

class Programmable(mixin.Named):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        includes_collective_effects: bool = False,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A programmable beam optics element.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    @property
    def beam_particles(
        self,
    ) -> collections.abc.Callable[
        [impactx.impactx_pybind.ImpactXParIter, impactx.impactx_pybind.RefPart], None
    ]:
        """
        hook for beam particles (pti, RefPart)
        """
    @beam_particles.setter
    def beam_particles(
        self,
        arg1: collections.abc.Callable[
            [impactx.impactx_pybind.ImpactXParIter, impactx.impactx_pybind.RefPart],
            None,
        ],
    ) -> None: ...
    @property
    def ds(self) -> float: ...
    @ds.setter
    def ds(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def includes_collective_effects(self) -> bool:
        """
        the push hooks already model collective effects over the length of this element: if True, no additional space charge, wakefield, CSR or ISR kicks are applied by ImpactX, default=False
        """
    @includes_collective_effects.setter
    def includes_collective_effects(self, arg1: bool) -> None: ...
    @property
    def nslice(self) -> int: ...
    @nslice.setter
    def nslice(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def push(
        self,
    ) -> collections.abc.Callable[
        [
            impactx.impactx_pybind.ImpactXParticleContainer,
            typing.SupportsInt | typing.SupportsIndex,
            typing.SupportsInt | typing.SupportsIndex,
        ],
        None,
    ]:
        """
        hook for push of whole container (pc, step, period)
        """
    @push.setter
    def push(
        self,
        arg1: collections.abc.Callable[
            [impactx.impactx_pybind.ImpactXParticleContainer, int, int], None
        ],
    ) -> None: ...
    @property
    def ref_particle(
        self,
    ) -> collections.abc.Callable[[impactx.impactx_pybind.RefPart], None]:
        """
        hook for reference particle (RefPart)
        """
    @ref_particle.setter
    def ref_particle(
        self, arg1: collections.abc.Callable[[impactx.impactx_pybind.RefPart], None]
    ) -> None: ...
    @property
    def threadsafe(self) -> bool:
        """
        allow threading via OpenMP for the particle iterator loop, default=False (note: if OMP backend is active)
        """
    @threadsafe.setter
    def threadsafe(self, arg1: bool) -> None: ...

class Quad(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        k: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        A Quadrupole magnet.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def k(self) -> float:
        """
        Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class RFCavity(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        *,
        ds,
        escale,
        freq,
        phase,
        cos_coefficients=None,
        sin_coefficients=None,
        z=None,
        field_on_axis=None,
        ncoef=None,
        dx=0,
        dy=0,
        rotation=0,
        aperture_x=0,
        aperture_y=0,
        mapsteps=10,
        nslice=1,
        name=None,
    ):
        """
        __init__(self: impactx.impactx_pybind.elements.RFCavity, ds: typing.SupportsFloat | typing.SupportsIndex, escale: typing.SupportsFloat | typing.SupportsIndex, freq: typing.SupportsFloat | typing.SupportsIndex, phase: typing.SupportsFloat | typing.SupportsIndex, cos_coefficients: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], sin_coefficients: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], dx: typing.SupportsFloat | typing.SupportsIndex = 0.0, dy: typing.SupportsFloat | typing.SupportsIndex = 0.0, rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0, aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0, aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0, mapsteps: typing.SupportsInt | typing.SupportsIndex = 10, nslice: typing.SupportsInt | typing.SupportsIndex = 1, name: str | None = None) -> None

        An RF cavity.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def set_coefficients(
        self,
        cos_coefficients: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
        sin_coefficients: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
    ) -> None:
        """
        Replace both arrays at once.

        They are validated together and applied as one step, so a rejected update
        leaves the element unchanged.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def cos_coefficients(self) -> list[float]:
        """
        cosine coefficients in the Fourier expansion of the on-axis RF electric field

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @cos_coefficients.setter
    def cos_coefficients(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def escale(self) -> float:
        """
        scaling factor for on-axis RF electric field in 1/m = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV)
        """
    @escale.setter
    def escale(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def freq(self) -> float:
        """
        RF frequency in Hz
        """
    @freq.setter
    def freq(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def map(self) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        linearized transport map around the reference particle (valid after a reference-particle push)
        """
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def phase(self) -> float:
        """
        RF driven phase in degrees
        """
    @phase.setter
    def phase(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def sin_coefficients(self) -> list[float]:
        """
        sine coefficients in the Fourier expansion of the on-axis RF electric field

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @sin_coefficients.setter
    def sin_coefficients(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def spin_coupling(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double:
        """
        linearized spin-orbit coupling matrix (valid after a reference-particle push)
        """

class Sbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        rc: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal sector bend.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def rc(self, ref: impactx.impactx_pybind.RefPart = None) -> float:
        """
        Radius of curvature in m
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """

class CFbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        rc: typing.SupportsFloat | typing.SupportsIndex,
        k: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal combined function bend (sector bend with quadrupole component).
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def k(self) -> float:
        """
        Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def rc(self) -> float:
        """
        Radius of curvature in m
        """
    @rc.setter
    def rc(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class Buncher(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        V: typing.SupportsFloat | typing.SupportsIndex,
        k: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A short linear RF cavity element at zero-crossing for bunching.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def V(self) -> float:
        """
        Normalized RF voltage drop V = Emax*L/(c*Brho)
        """
    @V.setter
    def V(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def k(self) -> float:
        """
        Wavenumber of RF in 1/m
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class ShortRF(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        V: typing.SupportsFloat | typing.SupportsIndex,
        freq: typing.SupportsFloat | typing.SupportsIndex,
        phase: typing.SupportsFloat | typing.SupportsIndex = -90.0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A short RF cavity element.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def V(self) -> float:
        """
        Normalized RF voltage V = maximum energy gain/(m*c^2)
        """
    @V.setter
    def V(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def freq(self) -> float:
        """
        RF frequency in Hz
        """
    @freq.setter
    def freq(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def phase(self) -> float:
        """
        RF synchronous phase in degrees (phase = 0 corresponds to maximum energy gain, phase = -90 corresponds go zero energy gain for bunching)
        """
    @phase.setter
    def phase(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class SoftSolenoid(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        *,
        ds,
        bscale,
        cos_coefficients=None,
        sin_coefficients=None,
        z=None,
        field_on_axis=None,
        ncoef=None,
        unit=0,
        dx=0,
        dy=0,
        rotation=0,
        aperture_x=0,
        aperture_y=0,
        mapsteps=10,
        nslice=1,
        name=None,
    ):
        """
        __init__(self: impactx.impactx_pybind.elements.SoftSolenoid, ds: typing.SupportsFloat | typing.SupportsIndex, bscale: typing.SupportsFloat | typing.SupportsIndex, cos_coefficients: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], sin_coefficients: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], unit: typing.SupportsInt | typing.SupportsIndex = 0, dx: typing.SupportsFloat | typing.SupportsIndex = 0.0, dy: typing.SupportsFloat | typing.SupportsIndex = 0.0, rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0, aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0, aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0, mapsteps: typing.SupportsInt | typing.SupportsIndex = 10, nslice: typing.SupportsInt | typing.SupportsIndex = 1, name: str | None = None) -> None

        A soft-edge solenoid.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def set_coefficients(
        self,
        cos_coefficients: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
        sin_coefficients: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
    ) -> None:
        """
        Replace both arrays at once.

        They are validated together and applied as one step, so a rejected update
        leaves the element unchanged.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def bscale(self) -> float:
        """
        Scaling factor for on-axis magnetic field Bz in inverse meters (if unit = 0) or magnetic field Bz in T (SI units, if unit = 1)
        """
    @bscale.setter
    def bscale(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def cos_coefficients(self) -> list[float]:
        """
        cosine coefficients in the Fourier expansion of the on-axis magnetic field Bz

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @cos_coefficients.setter
    def cos_coefficients(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def map(self) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        linearized transport map around the reference particle (valid after a reference-particle push)
        """
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def sin_coefficients(self) -> list[float]:
        """
        sine coefficients in the Fourier expansion of the on-axis magnetic field Bz

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @sin_coefficients.setter
    def sin_coefficients(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def spin_coupling(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double:
        """
        linearized spin-orbit coupling matrix (valid after a reference-particle push)
        """
    @property
    def spin_rotation_vector(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double:
        """
        reference spin rotation vector (valid after a reference-particle push)
        """
    @property
    def unit(self) -> int:
        """
        specification of units for scaling of the on-axis longitudinal magnetic field
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class Source(mixin.Named, mixin.Thin):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        distribution: str,
        openpmd_path: str,
        active_once: bool = True,
        load_ref_particle: bool = True,
        load_step: typing.SupportsInt | typing.SupportsIndex | None = None,
        load_step_index: typing.SupportsInt | typing.SupportsIndex | None = None,
        name: str | None = None,
    ) -> None:
        """
        A particle source.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def active_once(self) -> bool:
        """
        Inject particles only for the first lattice period.
        """
    @active_once.setter
    def active_once(self, arg1: bool) -> None: ...
    @property
    def distribution(self) -> str:
        """
        Distribution type of particles in the source
        """
    @distribution.setter
    def distribution(self, arg1: str) -> None: ...
    @property
    def load_ref_particle(self) -> bool:
        """
        Restore the reference particle from the species metadata of the openPMD file (particle tracking only).
        """
    @load_ref_particle.setter
    def load_ref_particle(self, arg1: bool) -> None: ...
    @property
    def load_step(self) -> int | None:
        """
        Which step (iteration) to load from the openPMD series: the ImpactX step at which the beam monitor wrote the beam, which is stored as the openPMD iteration in the file. Set at most one of load_step and load_step_index; if neither is set, the last step in the file is loaded.
        """
    @load_step.setter
    def load_step(
        self, arg1: typing.SupportsInt | typing.SupportsIndex | None
    ) -> None: ...
    @property
    def load_step_index(self) -> int | None:
        """
        Which step (iteration) to load from the openPMD series, by position in the file: 0 is the first step and -1 the last, counting back from it as in Python (-2 is the second to last step). Set at most one of load_step and load_step_index; if neither is set, the last step in the file is loaded.
        """
    @load_step_index.setter
    def load_step_index(
        self, arg1: typing.SupportsInt | typing.SupportsIndex | None
    ) -> None: ...
    @property
    def series_name(self) -> str:
        """
        Path to openPMD series as accepted by openPMD_api.Series
        """
    @series_name.setter
    def series_name(self, arg1: str) -> None: ...

class Sol(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        ds: typing.SupportsFloat | typing.SupportsIndex,
        ks: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        nslice: typing.SupportsInt | typing.SupportsIndex = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal hard-edge Solenoid magnet.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def ks(self) -> float:
        """
        Solenoid strength in m^(-1) (MADX convention) in (magnetic field Bz in T) / (rigidity in T-m)
        """
    @ks.setter
    def ks(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class PRot(mixin.Named, mixin.Thin):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        phi_in: typing.SupportsFloat | typing.SupportsIndex,
        phi_out: typing.SupportsFloat | typing.SupportsIndex,
        name: str | None = None,
    ) -> None:
        """
        An exact pole-face rotation in the x-z plane. Both angles are in degrees.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self, in_degrees: bool = False
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def phi_in(self) -> float:
        """
        angle of the reference particle with respect to the longitudinal (z) axis in the original frame in radian
        """
    @phi_in.setter
    def phi_in(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def phi_out(self) -> float:
        """
        angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in radian
        """
    @phi_out.setter
    def phi_out(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class SoftQuadrupole(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        *,
        ds,
        gscale,
        cos_coefficients=None,
        sin_coefficients=None,
        z=None,
        gradient_on_axis=None,
        ncoef=None,
        dx=0,
        dy=0,
        rotation=0,
        aperture_x=0,
        aperture_y=0,
        mapsteps=10,
        nslice=1,
        name=None,
    ):
        """
        __init__(self: impactx.impactx_pybind.elements.SoftQuadrupole, ds: typing.SupportsFloat | typing.SupportsIndex, gscale: typing.SupportsFloat | typing.SupportsIndex, cos_coefficients: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], sin_coefficients: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], dx: typing.SupportsFloat | typing.SupportsIndex = 0.0, dy: typing.SupportsFloat | typing.SupportsIndex = 0.0, rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0, aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0, aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0, mapsteps: typing.SupportsInt | typing.SupportsIndex = 10, nslice: typing.SupportsInt | typing.SupportsIndex = 1, name: str | None = None) -> None

        A soft-edge quadrupole.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def set_coefficients(
        self,
        cos_coefficients: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
        sin_coefficients: collections.abc.Sequence[
            typing.SupportsFloat | typing.SupportsIndex
        ],
    ) -> None:
        """
        Replace both arrays at once.

        They are validated together and applied as one step, so a rejected update
        leaves the element unchanged.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def cos_coefficients(self) -> list[float]:
        """
        cosine coefficients in the Fourier expansion of the on-axis field gradient

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @cos_coefficients.setter
    def cos_coefficients(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def gscale(self) -> float:
        """
        Scaling factor for on-axis field gradient in inverse meters
        """
    @gscale.setter
    def gscale(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def map(self) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        linearized transport map around the reference particle (valid after a reference-particle push)
        """
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def sin_coefficients(self) -> list[float]:
        """
        sine coefficients in the Fourier expansion of the on-axis field gradient

        Both arrays are always the same length. To change the length, set both at
        once with ``set_coefficients()``.
        """
    @sin_coefficients.setter
    def sin_coefficients(
        self,
        arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    ) -> None: ...
    @property
    def spin_coupling(
        self,
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double:
        """
        linearized spin-orbit coupling matrix (valid after a reference-particle push)
        """

class ThinDipole(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        theta: typing.SupportsFloat | typing.SupportsIndex,
        rc: typing.SupportsFloat | typing.SupportsIndex,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A thin kick model of a dipole bend.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self, in_degrees: bool = False
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def theta(self) -> float:
        """
        Bend angle (radian)
        """
    @theta.setter
    def theta(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class TaperedPL(mixin.Named, mixin.Thin, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        k: typing.SupportsFloat | typing.SupportsIndex,
        taper: typing.SupportsFloat | typing.SupportsIndex,
        unit: typing.SupportsInt | typing.SupportsIndex = 0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        A thin nonlinear plasma lens with transverse (horizontal) taper

                     .. math::

                        B_x = g \\left( y + \\frac{xy}{D_x} \\right), \\quad \\quad B_y = -g \\left(x + \\frac{x^2 + y^2}{2 D_x} \\right)

                     where :math:`g` is the (linear) field gradient in T/m and :math:`D_x` is the targeted horizontal dispersion in m.
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def k(self) -> float:
        """
        integrated focusing strength in m^(-1) (if unit = 0) or integrated focusing strength in T (if unit = 1)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def taper(self) -> float:
        """
        horizontal taper parameter in m^(-1) = 1 / (target horizontal dispersion in m)
        """
    @taper.setter
    def taper(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def unit(self) -> int:
        """
        specification of units for plasma lens focusing strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class LinearMap(mixin.Named, mixin.Alignment, mixin.PipeAperture):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        R: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ds: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_x: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        aperture_y: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        (A user-provided linear map, represented as a 6x6 transport matrix.)
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def R(self) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        linear map as a 6x6 transport matrix
        """
    @R.setter
    def R(
        self, arg1: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
    ) -> None: ...
    @property
    def ds(self) -> float:
        """
        segment length in m
        """
    @ds.setter
    def ds(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def nslice(self) -> int:
        """
        one, because we do not support slicing of this element
        """
    @property
    def symplectic(self) -> bool:
        """
        Check if the transport map is symplectic.

        A matrix R is symplectic if R^T J R = J, where J is the
        standard 6x6 skew-symmetric symplectic form (also called Omega).
        """

class SpinMap(mixin.Named, mixin.Alignment):
    def __eq__(self, other):
        """
        Value-based equality.

        Two elements are equal iff they are instances of the same class and
        their ``to_dict()`` outputs match key-for-key. Float values use
        ``==`` (so ``NaN != NaN``); list values are compared element-wise;
        AMReX SmallMatrix values use ``numpy.array_equal``.

        Returns ``NotImplemented`` when ``other`` is not the same element
        type, so Python's reflected-equality fallback applies (e.g., a
        foreign type's ``__eq__`` gets a chance, ultimately falling back to
        identity which yields ``False``).
        """
    def __hash__(self):
        """
        Value-based hash, consistent with ``__eq__``.

        Two elements that compare equal under ``==`` produce the same hash.
        The hash reflects the element's *current* parameter values; mutating
        an element after using it as a ``set`` member or ``dict`` key
        invalidates the container, the same contract a hashable mutable
        object would have.
        """
    def __init__(
        self,
        v: amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double,
        A: amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double,
        ds: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        rotation: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        name: str | None = None,
    ) -> None:
        """
        (A user-provided spin map, represented as a 3-vector and a 3x6 coupling matrix.)
        """
    def __repr__(self) -> str: ...
    def copy(self, **kwargs) -> typing.Any:
        """
        Return a new element with the same configuration.

        ``lattice.append(q)`` places ``q`` itself; ``lattice.append(q.copy())``
        places a separate element.

        Keyword arguments give the copy a different value for a parameter, which is
        how one element supplies the values for many::

            cavities = [rf.copy(name=f"rf{i}") for i in range(1, 5)]
            scan = [quad.copy(k=k) for k in (0.8, 0.9, 1.0)]

        Parameters that only mean something as a pair -- the cosine and sine
        coefficients, the normal and skew coefficients, the polygon vertices -- are
        given together, as the constructor takes them.
        """
    def finalize(self) -> None: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant equality for lattice elements.

        Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
        fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
        lists of floats are compared element-wise; AMReX SmallMatrix values
        use ``numpy.allclose``. All other value types (ints, strings,
        ``None``) fall back to strict ``==``. Mismatched element types and
        foreign operands return ``False`` (unless ``"type"`` is listed in
        ``ignore_attributes`` — see below).

        Parameters
        ----------
        other
            Element to compare against.
        rtol : float, optional
            Relative tolerance forwarded to ``math.isclose`` /
            ``numpy.allclose``. Default is ``1e-12`` — matches the
            ``dicts_equal`` helper used in the serialization tests, stricter
            than the ``math.isclose`` and ``numpy.isclose`` defaults.
        atol : float, optional
            Absolute tolerance. Default is ``0.0``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing. Useful when comparing
            loaded files where some bookkeeping fields (``"name"``) or even
            the element variant (``"type"``) should not affect the verdict.

            Including ``"type"`` disables the strict same-class check, so
            e.g. ``Drift`` and ``ExactDrift`` can be compared on their
            common parameters. Any remaining keys must still match between
            the two ``to_dict()`` outputs.
        """
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt | typing.SupportsIndex = 0,
        period: typing.SupportsInt | typing.SupportsIndex = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def reverse(self) -> None:
        """
        Reverse the element in-place so that pushing particles through
        it reverses the effect of the original element.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | bool
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
        | None,
    ]: ...
    def transfer_map(
        self, ref: impactx.impactx_pybind.RefPart
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Return this element's 6x6 linear transport map for the given
        reference particle.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).
        For an element with ``nslice`` > 1 this is the map of a single
        ``ds/nslice`` slice (the building block that the lattice transfer
        map composes). Raises for an element whose linear transport map is
        not implemented.

        Uses the same convention as ``KnownElementsList.transfer_map``: the
        reference particle is advanced through one slice first, and the map
        is evaluated at that post-advance reference state. The reference
        particle you pass in is copied and left unmodified.

        :param ref: reference particle at the element entrance
        """
    @property
    def A(self) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double:
        """
        spin-orbit coupling generator of rotation as a 3x6 matrix
        """
    @A.setter
    def A(
        self, arg1: amrex.space3d.amrex_3d_pybind.SmallMatrix_3x6_F_SI1_double
    ) -> None: ...
    @property
    def ds(self) -> float:
        """
        segment length in m
        """
    @ds.setter
    def ds(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def nslice(self) -> int:
        """
        one, because we do not support slicing of this element
        """
    @property
    def v(self) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double:
        """
        design axis-angle generator of spin rotation as a 3x1 vector
        """
    @v.setter
    def v(
        self, arg1: amrex.space3d.amrex_3d_pybind.SmallMatrix_3x1_F_SI1_double
    ) -> None: ...

class KnownElementsList:
    def __contains__(self, element: typing.Any) -> bool: ...
    @typing.overload
    def __delitem__(self, slice: slice) -> None:
        """
        Remove the selected elements.
        """
    @typing.overload
    def __delitem__(self, index: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Remove the element at a position.

        An element that also sits at other positions keeps those.
        """
    def __eq__(self, other):
        """
        Element-wise equality with any iterable of elements.

        Duck-typed across container types: a ``KnownElementsList`` compares
        equal to a plain Python ``list`` of elements or to a
        ``FilteredElementsList`` as long as their lengths match and every
        pair of elements compares equal under ``==``. Returns
        ``NotImplemented`` for non-iterable operands so Python's
        reflected-equality fallback applies. Mutable containers are
        deliberately unhashable (``__hash__ = None``), matching the Python
        ``list`` convention.
        """
    @typing.overload
    def __getitem__(
        self, index: typing.SupportsInt | typing.SupportsIndex
    ) -> typing.Any:
        """
        Return the element at a position.
        """
    @typing.overload
    def __getitem__(self, slice: slice) -> typing.Any:
        """
        Return a new lattice holding the selected elements.
        """
    def __init__(self, elements=None):
        """
        Create a lattice, optionally filled with elements.

        ``elements`` is a single element or any iterable of elements: a list, another
        lattice, a selection, a generator. The elements are shared, not copied, so the
        caller keeps handles to the very elements the lattice holds. Constructing goes
        through ``extend`` for exactly that reason: the C++ constructor cannot see the
        object being constructed, and so cannot record which Python objects own the
        elements.
        """
    def __iter__(self) -> typing.Any:
        """
        Iterate over the elements of the lattice.
        """
    def __len__(self) -> int:
        """
        The length of the list.
        """
    def __reversed__(self) -> typing.Any:
        """
        Iterate over the elements from the end.
        """
    @typing.overload
    def __setitem__(
        self, index: typing.SupportsInt | typing.SupportsIndex, element: typing.Any
    ) -> None:
        """
        Replace the element at a position.
        """
    @typing.overload
    def __setitem__(self, slice: slice, elements: collections.abc.Iterable) -> None:
        """
        Replace the selected elements.
        """
    def append(self, element: typing.Any) -> None:
        """
        Add an element to the end of the lattice.

        The lattice holds this element, so changing it afterwards changes what is
        tracked. Adding the same element twice places it at two positions. For two
        independent elements, add ``element.copy()`` or construct a second one.
        """
    def clear(self) -> None:
        """
        Remove all elements from the lattice.
        """
    def count(self, element: typing.Any) -> int:
        """
        Return how many positions this element occupies.
        """
    def count_by_kind(self, kind_pattern) -> int:
        """
        Count elements of a specific kind.

        Args:
            kind_pattern: The element kind to count. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            int: Number of elements of the specified kind.
        """
    def extend(self, elements: collections.abc.Iterable) -> typing.Any:
        """
        Add several elements to the lattice.
        """
    def from_dicts(self, dicts: list[dict]):
        """
        Load and append elements from a list of dictionaries.

        Each dictionary should be in the format produced by element.to_dict(),
        containing at minimum a 'type' key identifying the element class.

        Args:
            dicts: List of element dictionaries

        Example:
            .. code-block:: python

                import json
                from impactx import elements

                # Load from JSON
                with open("lattice.impactx.json") as f:
                    data = json.load(f)

                lattice = elements.KnownElementsList()
                lattice.from_dicts(data)

        Note:
            Elements with matrix parameters (LinearMap, SpinMap) require
            the matrices to be AMReX SmallMatrix objects. Use
            :func:`impactx.extensions.matrix_hook` as a JSON object_hook
            when loading such elements.
        """
    def from_pals(self, pals_beamline, nslice=1, *, min_model="linear"):
        """
        Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) object.

        ``min_model`` is the element model floor, see :py:func:`load_file`.

        https://github.com/campa-consortium/pals-python
        """
    def get_kinds(self) -> list[type]:
        """
        Get all unique element kinds in the list.

        Returns:
            list[type]: List of unique element types (sorted by name).
        """
    def has_kind(self, kind_pattern) -> bool:
        """
        Check if list contains elements of a specific kind.

        Args:
            kind_pattern: The element kind to check for. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            bool: True if at least one element of the specified kind exists.
        """
    def index(self, element: typing.Any) -> int:
        """
        Return the first position holding this element.
        """
    def insert(
        self, index: typing.SupportsInt | typing.SupportsIndex, element: typing.Any
    ) -> None:
        """
        Insert an element before a position.
        """
    def is_empty(self) -> bool: ...
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant element-wise comparison with any iterable of elements.

        For each pair of elements at the same index, calls the element's own
        ``isclose(rtol=..., atol=..., ignore_attributes=...)``. Lengths must
        match; foreign or non-iterable operands return ``False``. Defaults
        match the per-element ``isclose`` (``rtol=1e-12``, ``atol=0.0``).

        Parameters
        ----------
        other : iterable of elements
            Any container (``KnownElementsList``, ``FilteredElementsList``,
            or plain ``list``) whose elements expose ``isclose``.
        rtol, atol : float
            Forwarded to each element's ``isclose``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing each pair of elements.
            Includes the special key ``"type"`` to compare across element
            variants (e.g., ``Drift`` vs ``ExactDrift``). Forwarded to each
            element's ``isclose``.
        """
    def load_file(self, filename, nslice=1, *, min_model="linear"):
        """
        Load and append a lattice file from MAD-X (.madx) or PALS (e.g., .pals.yaml) formats.

        ``min_model`` raises the lower gate of the element model selection: the reader
        still picks the cheapest ImpactX element that represents the imported element,
        but never one below the requested tier (``"linear"``, ``"paraxial"`` or
        ``"exact"``). Where a tier is not implemented for an element family, the next
        higher one is used. Where no model reaches the requested tier at all, as for a
        solenoid, which has no exact model, a warning is emitted.

        .. warning::

           Our MAD-X and PALS parsers are under active development
           and provided as a preview. Please check any loaded lattice
           files very carefully. Please report your experience and bugs
           on our `issue tracker <https://github.com/BLAST-ImpactX/impactx/issues>`__.
        """
    def map_trace(self, ref: impactx.impactx_pybind.RefPart) -> list:
        """
        Trace the cumulative 6x6 linear transport map element by element.

        The reference particle is passed by value (intentional copy); the
        caller's reference particle is not modified in place. This matches
        the convention used by ``transfer_map``.

        This per-element trace is what ``sim.twiss()`` consumes to transport
        Twiss functions through the lattice.

        If you only need the final cumulative map at the lattice exit,
        prefer ``transfer_map(ref)`` instead of indexing the last entry
        of ``map_trace(ref)``.

        :param ref: A reference particle.
        :return: A list of dictionaries, one per lattice element plus a
            leading entry for the starting position. Each entry contains:

            * ``s``    -- integrated path length along the reference orbit,
              in meters;
            * ``name`` -- user-supplied element name (empty string if not
              named);
            * ``type`` -- element type string (e.g. ``"Drift"``,
              ``"Quad"``, ``"Sbend"``);
            * ``M``    -- cumulative 6x6 linear transport map from the
              start of the lattice to the exit of this element (a
              ``Map6x6`` instance; call ``.to_numpy()`` for a standard
              C-ordered NumPy array).

            The first entry always has the identity map at the starting
            ``s``; the last entry contains the same map as ``transfer_map``.
        """
    def plot_survey(
        self, ref=None, ax=None, legend=True, legend_ncols=5, palette="cern-lhc"
    ):
        """
        Plot over s of all elements in the KnownElementsList.

        A positive element strength denotes horizontal focusing (e.g. for quadrupoles) and bending to the right (for dipoles).  In general, this depends on both the sign of the field and the sign of the charge.

        Parameters
        ----------
        self : ImpactXParticleContainer_*
            The KnownElementsList class in ImpactX
        ref : RefPart
            A reference particle, checked for the charge sign to plot focusing/defocusing strength directions properly.
        ax : matplotlib axes
            A plotting area in matplotlib (called axes there).
        legend: bool
            Plot a legend if true.
        legend_ncols: int
            Number of columns for lattice element types in the legend.
        palette: string
            Color palette.

        Returns
        -------
        Either populates the matplotlib axes in ax or creates a new axes containing the plot.
        """
    def pop_back(self) -> typing.Any:
        """
        Remove and return the last element of the lattice.
        """
    def remove(self, element: typing.Any) -> None:
        """
        Remove the first occurrence of this element.
        """
    def select(self, *, kind=None, name=None) -> FilteredElementsList:
        """
        Filter elements by type and name with OR-based logic.

        This method supports filtering elements by their type and/or name using keyword arguments.
        Returns references to original elements, allowing modification and chaining.

        **Filtering Logic:**

        - **Within a single filter**: OR logic (e.g., ``kind=["Drift", "Quad"]`` matches Drift OR Quad)
        - **Between different filters**: OR logic (e.g., ``kind="Quad", name="quad1"`` matches Quad OR named "quad1")
        - **Chaining filters**: AND logic (e.g., ``lattice.select(kind="Drift").select(name="drift1")`` matches Drift AND named "drift1")

        :param kind: Element type(s) to filter by. Can be a single string/type or a list/tuple
                     of strings/types for OR-based filtering. String values support exact matches
                     and regex patterns. Examples: "Drift", r".*Quad", elements.Drift, ["Drift", r".*Quad"], [elements.Drift, elements.Quad]
        :type kind: str or type or list[str | type] or tuple[str | type, ...] or None, optional

        :param name: Element name(s) to filter by. Can be a single string, regex pattern string, or
                     a list/tuple of strings and/or regex pattern strings for OR-based filtering.
                     Examples: "quad1", r"quad\\d+", ["quad1", "quad2"], [r"quad\\d+", "bend1"]
        :type name: str or list[str] or tuple[str, ...] or None, optional

        :return: FilteredElementsList containing references to original elements
        :rtype: FilteredElementsList

        :raises TypeError: If kind/name parameters have wrong types

        **Examples:**

        Single value filtering:

        .. code-block:: python

            lattice.select(kind="Drift")  # Get all drift elements (string)
            lattice.select(kind=elements.Drift)  # Get all drift elements (type)
            lattice.select(
                kind=r".*Quad"
            )  # Get all elements matching regex pattern (Quad, ExactQuad, ChrQuad)
            lattice.select(name="quad1")  # Get elements named "quad1"
            lattice.select(
                kind="Quad", name="quad1"
            )  # Get quad elements OR elements named "quad1"

        OR-based filtering with lists (within single filter):

        .. code-block:: python

            lattice.select(kind=["Drift", "Quad"])  # Get drift OR quad elements (strings)
            lattice.select(kind=[elements.Drift, elements.Quad])  # Get drift OR quad elements (types)
            lattice.select(kind=["Drift", elements.Quad])  # Mix strings and types
            lattice.select(kind=[r".*Quad", r".*Bend.*"])  # Mix regex patterns
            lattice.select(name=["quad1", "quad2"])  # Get elements named "quad1" OR "quad2"

         Regex pattern filtering:

         .. code-block:: python

             lattice.select(name=r"quad\\d+")  # Get elements matching pattern
             lattice.select(name=[r"quad\\d+", "bend1"])  # Mix regex and strings

        Chaining filters (AND logic between chained calls):

        .. code-block:: python

            lattice.select(kind="Drift").select(
                name="drift1"
            )  # Drift elements AND named "drift1"
            lattice.select(kind="Quad")[0]  # First quad element
            lattice.select(name="quad1").select(
                kind="Quad"
            )  # Elements named "quad1" AND of type "Quad"

        Reference preservation and modification:

        .. code-block:: python

            drift_elements = lattice.select(kind="Drift")
            drift_elements[0].ds = 5.0  # Modifies the original element in lattice
            assert lattice[0].ds == 5.0  # Original element is modified

        Modification of elements (reference preservation):

        .. code-block:: python

            drift = lattice.select(kind="Drift")[0]  # Get first drift element
            drift.ds = 2.0  # Modify original element
            quad_elements = lattice.select(kind="Quad")  # Get all quad elements
            quad_elements[0].k = 1.5  # Modify first quad's strength
            # All modifications affect the original lattice elements
        """
    def size(self) -> int: ...
    def to_dicts(self) -> list[dict]:
        """
        Serialize the lattice to a list of dictionaries.

        Each element is converted to a dictionary using its to_dict() method.
        The resulting list can be serialized to JSON, YAML, or other formats.

        Returns:
            list[dict]: List of element dictionaries

        Example:
            .. code-block:: python

                import json
                from impactx import elements

                lattice = elements.KnownElementsList(
                    [
                        elements.Drift(ds=1.0, name="d1"),
                        elements.Quad(ds=0.5, k=2.0, name="q1"),
                    ]
                )

                # Serialize to JSON
                data = lattice.to_dicts()
                with open("lattice.impactx.json", "w") as f:
                    json.dump(data, f, indent=2)

        Note:
            Elements with matrix parameters (LinearMap, SpinMap) contain
            AMReX SmallMatrix objects that require custom JSON encoding.
            Use :func:`impactx.extensions.ImpactXEncoder` for JSON serialization
            of such elements.
        """
    def to_py(self) -> str:
        """
        Generate Python code that recreates this lattice.

        Returns a string containing a complete Python script with imports
        and a ``get_lattice()`` function that returns a KnownElementsList
        with all elements.

        Returns:
            str: Python source code

        Example:
            .. code-block:: python

                from impactx import elements

                lattice = elements.KnownElementsList(
                    [
                        elements.Drift(ds=1.0, name="d1"),
                        elements.Quad(ds=0.5, k=2.0, name="q1"),
                    ]
                )

                # Generate Python code
                code = lattice.to_py()
                print(code)

                # Save to file
                with open("my_lattice.py", "w") as f:
                    f.write(code)

                # Later, use the generated file:
                # from my_lattice import get_lattice
                # lattice = get_lattice()
        """
    def transfer_map(
        self,
        ref: impactx.impactx_pybind.RefPart,
        order: str = "linear",
        fallback_identity_map: bool = False,
    ) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        Calculate the end-to-end transfer map of the elements in the list.

        Currently only the linear transfer map is implemented (``order="linear"``);
        the ``order`` parameter is reserved for future higher-order extensions.
        In linear mode the 6x6 map is composed element by element, using each
        element's analytic per-slice linear transport map.

        Collective effects like space charge, Coherent/Incoherent Synchrotron
        Radiation (CSR/ISR), and wakefield effects are not applied here; the
        returned map describes the purely linear single-particle dynamics of the
        design lattice.

        Phase-space ordering in the returned matrix is (x, px, y, py, t, pt).

        :param ref: reference particle at the starting s
        :param order: So far, only the calculation of linear transfer maps is supported.
        :param fallback_identity_map: For elements with an undefined transfer map in the lattice, assume the identity matrix.
        """
    @property
    def generation(self) -> int:
        """
        How often the sequence of elements changed.

        Counts structural edits only: changing a parameter on an element that is
        already in the lattice does not move anything and does not change this.
        A view that remembers positions compares this to notice they went stale.

        Compare it for equality, not by how far it moved: one call can edit the
        sequence in several steps, as deleting a selection does.
        """

class FilteredElementsList:
    """
    Result of ``KnownElementsList.select(...)`` or chained ``.select()`` calls: a filtered
    view of the same underlying lattice.

    Indexing (``self[i]``) returns elements from the original ``KnownElementsList``; changing
    fields on those elements modifies the lattice in place. You can narrow the filter again with
    ``.select(...)`` (AND logic between chained calls). After ``delete``, ``replace_each``, or
    ``replace_with_drifts``, obtain a new selection from the lattice; earlier filter objects must
    not be used.
    """
    def __eq__(self, other):
        """
        Element-wise equality with any iterable of elements.

        Duck-typed across container types: a ``KnownElementsList`` compares
        equal to a plain Python ``list`` of elements or to a
        ``FilteredElementsList`` as long as their lengths match and every
        pair of elements compares equal under ``==``. Returns
        ``NotImplemented`` for non-iterable operands so Python's
        reflected-equality fallback applies. Mutable containers are
        deliberately unhashable (``__hash__ = None``), matching the Python
        ``list`` convention.
        """
    def __getitem__(self, key): ...
    def __init__(self, original_list, indices): ...
    def __iter__(self): ...
    def __len__(self): ...
    def __repr__(self): ...
    def __str__(self): ...
    def _require_valid(self) -> None:
        """
        Raise if the lattice changed after this view was taken.
        """
    def count_by_kind(self, kind_pattern) -> int:
        """
        Count elements of a specific kind in the filtered list.

        Args:
            kind_pattern: The element kind to count. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            int: Number of elements of the specified kind.
        """
    def delete(self) -> None:
        """
        Remove selected elements from the underlying lattice. Invalidates this and all other
        live selections on the same lattice. Returns None.
        """
    def get_kinds(self) -> list[type]:
        """
        Get all unique element kinds in the filtered list.

        Returns:
            list[type]: List of unique element types (sorted by name).
        """
    def has_kind(self, kind_pattern) -> bool:
        """
        Check if filtered list contains elements of a specific kind.

        Args:
            kind_pattern: The element kind to check for. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            bool: True if at least one element of the specified kind exists.
        """
    def isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
        """
        Tolerant element-wise comparison with any iterable of elements.

        For each pair of elements at the same index, calls the element's own
        ``isclose(rtol=..., atol=..., ignore_attributes=...)``. Lengths must
        match; foreign or non-iterable operands return ``False``. Defaults
        match the per-element ``isclose`` (``rtol=1e-12``, ``atol=0.0``).

        Parameters
        ----------
        other : iterable of elements
            Any container (``KnownElementsList``, ``FilteredElementsList``,
            or plain ``list``) whose elements expose ``isclose``.
        rtol, atol : float
            Forwarded to each element's ``isclose``.
        ignore_attributes : str or iterable of str, optional
            ``to_dict()`` keys to skip when comparing each pair of elements.
            Includes the special key ``"type"`` to compare across element
            variants (e.g., ``Drift`` vs ``ExactDrift``). Forwarded to each
            element's ``isclose``.
        """
    def replace_each(self, element, *, keep_name=True, keep_ds=False):
        """
        Replace each selected element with a copy of ``element``, optionally keeping name and
        ``ds`` from the replaced element (``keep_ds`` defaults to False). Invalidates prior views;
        returns a new selection over the same indices.
        """
    def replace_with_drifts(
        self, *, model="match", keep_alignment=True, keep_aperture=False
    ):
        """
        Replace each selected element with a drift of the matching physics family.

        When ``model="match"``: ``Exact*`` elements become ``ExactDrift``, ``Chr*`` elements
        become ``ChrDrift``, and all other (linear) elements become ``Drift``. When
        ``model`` is ``"linear"``, ``"paraxial"``, or ``"exact"``, every selected slot uses
        that drift model. Names and segment length ``ds`` are always taken from the replaced
        element.

        By default, alignment errors (dx, dy, rotation) are preserved and apertures are
        cleared. Use ``keep_alignment=False`` to zero alignment errors, or
        ``keep_aperture=True`` to preserve aperture_x/aperture_y.
        """
    def select(self, *, kind=None, name=None):
        """
        Apply filtering to this filtered list.

        This method applies additional filtering to an already filtered list,
        maintaining references to the original elements and enabling chaining.

        **Filtering Logic:**

        - **Within a single filter**: OR logic (e.g., ``kind=["Drift", "Quad"]`` matches Drift OR Quad)
        - **Between different filters**: OR logic (e.g., ``kind="Quad", name="quad1"`` matches Quad OR named "quad1")
        - **Chaining filters**: AND logic (e.g., ``lattice.select(kind="Drift").select(name="drift1")`` matches Drift AND named "drift1")

        :param kind: Element type(s) to filter by. Can be a single string/type or a list/tuple
                     of strings/types for OR-based filtering. String values support exact matches
                     and regex patterns. Examples: "Drift", r".*Quad", elements.Drift, ["Drift", r".*Quad"], [elements.Drift, elements.Quad]
        :type kind: str or type or list[str | type] or tuple[str | type, ...] or None, optional

        :param name: Element name(s) to filter by. Can be a single string, regex pattern string, or
                     a list/tuple of strings and/or regex pattern strings for OR-based filtering.
                     Examples: "quad1", r"quad\\d+", ["quad1", "quad2"], [r"quad\\d+", "bend1"]
        :type name: str or list[str] or tuple[str, ...] or None, optional

        :return: FilteredElementsList containing references to original elements
        :rtype: FilteredElementsList

        :raises TypeError: If kind/name parameters have wrong types

        **Examples:**

        Additional filtering on already filtered results:

        .. code-block:: python

            drift_elements = lattice.select(
                kind="Drift"
            )  # or lattice.select(kind=elements.Drift)
            first_drift = drift_elements.select(
                name="drift1"
            )  # Further filter drifts by name
            quad_elements = lattice.select(
                kind="Quad"
            )  # or lattice.select(kind=elements.Quad)
            strong_quads = quad_elements.select(
                name=r"quad\\d+"
            )  # Filter quads by regex pattern
        """
