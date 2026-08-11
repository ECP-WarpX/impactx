"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell, Edoardo Zoni
License: BSD-3-Clause-LBNL
"""

import os
import re

from impactx import Config, elements

from ..element_models import DRIFT_MODEL_CLASSES, tier_of_class, validate_model

FILTERED_ELEMENTS_LIST_INVALID_MSG = (
    "This lattice selection is no longer valid because the lattice was modified; "
    "call select() again on the lattice."
)


def _drift_class_for_replace_with_drifts(model: str, old_el) -> type:
    """Map ``model`` and ``old_el`` to the Drift / ChrDrift / ExactDrift class to insert.

    For ``model=="match"``, the class follows the replaced element's own tier (see
    ``impactx.element_models.tier_of_class``). Otherwise ``model`` must already be
    validated against :data:`impactx.element_models.MODEL_TIERS`."""
    if model == "match":
        key = tier_of_class(type(old_el).__name__)
    else:
        key = model
    return DRIFT_MODEL_CLASSES[key]


# These report an angle in radians from ``to_dict()`` although their constructors
# take degrees, so serializing them for round-trip needs ``in_degrees=True``.
# See https://github.com/BLAST-ImpactX/impactx/issues/1367
_DEGREE_ELEMENTS = (
    "ExactSbend",
    "PlaneXYRot",
    "PRot",
    "ThinDipole",
)


def _element_to_dict(element) -> dict:
    """Serialize an element to a dict that its constructor can consume again.

    Applies the radians/degrees work-around for the element types whose
    ``to_dict()`` disagrees with their constructor.
    """
    if type(element).__name__ in _DEGREE_ELEMENTS:
        return element.to_dict(in_degrees=True)
    return element.to_dict()


def _make_drift_from_old(
    cls,
    old_el,
    *,
    keep_name,
    keep_ds,
    keep_alignment,
    keep_aperture,
):
    """Build a drift (``cls`` is Drift / ChrDrift / ExactDrift) from thick-element fields on ``old_el``.

    When ``keep_ds`` is False, ``ds`` is set to 0. When ``keep_name`` is False, ``name`` is None.
    If ``keep_ds`` is True but ``old_el`` has no ``ds`` attribute (thin element), ``ds`` defaults to 0.
    When ``keep_alignment`` is True, copy dx/dy/rotation from the old element.
    When ``keep_aperture`` is True, copy aperture_x/aperture_y from the old element.
    """
    nslice = getattr(old_el, "nslice", 1)

    if keep_ds and hasattr(old_el, "ds"):
        ds = old_el.ds
    else:
        ds = 0.0

    name = None
    if keep_name:
        if hasattr(old_el, "has_name") and old_el.has_name:
            name = old_el.name

    if keep_alignment:
        dx = getattr(old_el, "dx", 0.0)
        dy = getattr(old_el, "dy", 0.0)
        rotation = getattr(old_el, "rotation", 0.0)
    else:
        dx = 0.0
        dy = 0.0
        rotation = 0.0

    if keep_aperture:
        aperture_x = getattr(old_el, "aperture_x", 0.0)
        aperture_y = getattr(old_el, "aperture_y", 0.0)
    else:
        aperture_x = 0.0
        aperture_y = 0.0

    return cls(
        ds=ds,
        dx=dx,
        dy=dy,
        rotation=rotation,
        aperture_x=aperture_x,
        aperture_y=aperture_y,
        nslice=nslice,
        name=name,
    )


def load_file(self, filename, nslice=1, *, min_model="linear"):
    """Load and append a lattice file from MAD-X (.madx) or PALS (e.g., .pals.yaml) formats.

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

    import warnings

    warnings.warn(
        "Our MAD-X and PALS parsers are under active development and provided as a "
        "preview. Please check any loaded lattice files very carefully. Please "
        "report your experience and bugs on our issue tracker: "
        "https://github.com/BLAST-ImpactX/impactx/issues",
        RuntimeWarning,
        stacklevel=2,
    )

    # Attempt to strip two levels of file extensions to determine the schema.
    #   Examples: fodo.madx, fodo.pals.yaml, fodo.pals.json, ...
    file_noext, extension = os.path.splitext(filename)
    file_noext_noext, extension_inner = os.path.splitext(file_noext)

    if extension == ".madx":
        # example: fodo.madx
        from ..madx_to_impactx import read_lattice

        # TODO: Expose explicit MAD-X line/sequence selection in this public API
        # once the user-facing interface is settled. The lower-level translator
        # already supports read_lattice(..., line=..., sequence=...).
        self.extend(
            read_lattice(
                filename, nslice, line=None, sequence=None, min_model=min_model
            )
        )
        return

    elif extension_inner == ".pals":
        from pals import load as load_pals_file

        self.from_pals(load_pals_file(filename), nslice, min_model=min_model)
        return

    raise RuntimeError(
        f"load_file: No support for file {filename} with extension {extension} yet."
    )


def from_pals(self, pals_beamline, nslice=1, *, min_model="linear"):
    """Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) object.

    ``min_model`` is the element model floor, see :py:func:`load_file`.

    https://github.com/campa-consortium/pals-python
    """
    from ..pals_to_impactx import read_lattice

    self.extend(read_lattice(pals_beamline, nslice, min_model=min_model))


class FilteredElementsList:
    """Result of ``KnownElementsList.select(...)`` or chained ``.select()`` calls: a filtered
    view of the same underlying lattice.

    Indexing (``self[i]``) returns elements from the original ``KnownElementsList``; changing
    fields on those elements modifies the lattice in place. You can narrow the filter again with
    ``.select(...)`` (AND logic between chained calls). After ``delete``, ``replace_each``, or
    ``replace_with_drifts``, obtain a new selection from the lattice; earlier filter objects must
    not be used.
    """

    def __init__(self, original_list, indices):
        self._original_list = original_list
        self._indices = list(indices) if not isinstance(indices, list) else indices
        # A selection is a list of positions, so any edit that moves elements makes it
        # describe something else. The lattice counts those edits, which catches every
        # edit, whether it was made on the lattice or through another selection.
        self._generation = original_list.generation

    def _require_valid(self) -> None:
        """Raise if the lattice changed after this view was taken."""
        if self._original_list.generation != self._generation:
            raise RuntimeError(FILTERED_ELEMENTS_LIST_INVALID_MSG)

    def __getitem__(self, key):
        self._require_valid()
        if isinstance(key, int):
            return self._original_list[self._indices[key]]
        elif isinstance(key, slice):
            sliced_indices = self._indices[key]
            return FilteredElementsList(self._original_list, sliced_indices)
        else:
            raise TypeError(f"Invalid key type: {type(key)}")

    def __len__(self):
        self._require_valid()
        return len(self._indices)

    def __iter__(self):
        self._require_valid()
        for i in self._indices:
            yield self._original_list[i]

    def select(
        self,
        *,
        kind=None,
        name=None,
    ):
        r"""Apply filtering to this filtered list.

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
                     Examples: "quad1", r"quad\d+", ["quad1", "quad2"], [r"quad\d+", "bend1"]
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
                name=r"quad\d+"
            )  # Filter quads by regex pattern
        """
        self._require_valid()
        # Apply filtering directly to the indices we already have
        if kind is not None or name is not None:
            # Validate parameters
            _validate_select_parameters(kind, name)

            matching_indices = []

            for i in self._indices:
                element = self._original_list[i]
                if _check_element_match(element, kind, name):
                    matching_indices.append(i)

            return FilteredElementsList(self._original_list, matching_indices)

        # If no filtering criteria provided, return all current elements
        return FilteredElementsList(self._original_list, self._indices)

    def delete(self) -> None:
        """Remove selected elements from the underlying lattice. Invalidates this and all other
        live selections on the same lattice. Returns None."""
        self._require_valid()
        original = self._original_list
        to_remove = set(self._indices)
        if not to_remove:
            # nothing selected changes nothing, so other selections stay usable
            return None

        # Keep the rest, in one pass. Deleting the selected positions one at a time moves
        # everything after each of them, which is quadratic in the length of the lattice.
        # Elements that were not selected are carried over as they are: they keep their
        # identity, their Python subclass and their callbacks.
        kept = [original[i] for i in range(len(original)) if i not in to_remove]
        original.clear()
        original.extend(kept)
        return None

    def replace_each(self, element, *, keep_name=True, keep_ds=False):
        """Replace each selected element with a copy of ``element``, optionally keeping name and
        ``ds`` from the replaced element (``keep_ds`` defaults to False). Invalidates prior views;
        returns a new selection over the same indices."""
        self._require_valid()
        original = self._original_list
        indices = self._indices
        if not indices:
            # nothing selected changes nothing, so other selections stay usable
            return FilteredElementsList(original, [])

        # One copy per selected position, so the replacements are independent elements
        # rather than one element repeated. Unselected positions are left alone.
        #
        # Build every replacement before installing any of them: a template that cannot be
        # copied, or whose name or length cannot be set, must leave the lattice as it was
        # rather than replacing the positions reached before it failed.
        replacements = []
        for i in indices:
            old_el = original[i]
            new_el = element.copy()
            if keep_name:
                if hasattr(old_el, "has_name") and old_el.has_name:
                    new_el.name = old_el.name
            if keep_ds and hasattr(old_el, "ds"):
                new_el.ds = old_el.ds
            replacements.append((i, new_el))

        for i, new_el in replacements:
            original[i] = new_el

        return FilteredElementsList(original, list(indices))

    def replace_with_drifts(
        self, *, model="match", keep_alignment=True, keep_aperture=False
    ):
        """Replace each selected element with a drift of the matching physics family.

        When ``model="match"``: ``Exact*`` elements become ``ExactDrift``, ``Chr*`` elements
        become ``ChrDrift``, and all other (linear) elements become ``Drift``. When
        ``model`` is ``"linear"``, ``"paraxial"``, or ``"exact"``, every selected slot uses
        that drift model. Names and segment length ``ds`` are always taken from the replaced
        element.

        By default, alignment errors (dx, dy, rotation) are preserved and apertures are
        cleared. Use ``keep_alignment=False`` to zero alignment errors, or
        ``keep_aperture=True`` to preserve aperture_x/aperture_y."""
        self._require_valid()
        original = self._original_list
        indices = self._indices
        if not indices:
            # nothing selected changes nothing, so other selections stay usable
            return FilteredElementsList(original, [])

        validate_model(model, argument="model", extra_values=("match",))

        # Only the selected positions are rewritten; everything else stays as it is. Every
        # drift is built before any is installed, so a failure part-way leaves the lattice
        # as it was.
        replacements = [
            (
                i,
                _make_drift_from_old(
                    _drift_class_for_replace_with_drifts(model, original[i]),
                    original[i],
                    keep_name=True,
                    keep_ds=True,
                    keep_alignment=keep_alignment,
                    keep_aperture=keep_aperture,
                ),
            )
            for i in indices
        ]

        for i, new_el in replacements:
            original[i] = new_el

        return FilteredElementsList(original, list(indices))

    def get_kinds(self) -> list[type]:
        """Get all unique element kinds in the filtered list.

        Returns:
            list[type]: List of unique element types (sorted by name).
        """
        self._require_valid()
        return get_kinds(self)

    def count_by_kind(self, kind_pattern) -> int:
        """Count elements of a specific kind in the filtered list.

        Args:
            kind_pattern: The element kind to count. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            int: Number of elements of the specified kind.
        """
        self._require_valid()
        return count_by_kind(self, kind_pattern)

    def has_kind(self, kind_pattern) -> bool:
        """Check if filtered list contains elements of a specific kind.

        Args:
            kind_pattern: The element kind to check for. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            bool: True if at least one element of the specified kind exists.
        """
        self._require_valid()
        return has_kind(self, kind_pattern)

    def __repr__(self):
        self._require_valid()
        return f"FilteredElementsList({len(self)} elements)"

    def __str__(self):
        self._require_valid()
        return f"FilteredElementsList({len(self)} elements)"


def _is_regex_pattern(pattern: str) -> bool:
    """Check if a string looks like a regex pattern by testing if it contains regex metacharacters."""
    # Simple heuristic: if it contains regex metacharacters, treat as regex
    regex_chars = r".*+?^${}[]|()\\"
    return any(char in pattern for char in regex_chars)


def _matches_string(text: str, string_pattern: str) -> bool:
    """Check if text matches a string pattern (exact match or regex)."""
    if _is_regex_pattern(string_pattern):
        try:
            return bool(re.search(string_pattern, text))
        except re.error:
            # If regex compilation fails, fall back to exact match
            return text == string_pattern
    else:
        return text == string_pattern


def _validate_select_parameters(kind, name):
    """Validate parameters for select methods.

    Args:
        kind: Element type(s) to filter by
        name: Element name(s) to filter by

    Raises:
        TypeError: If parameters have wrong types
    """
    if kind is not None:
        if isinstance(kind, (list, tuple)):
            for k in kind:
                if not isinstance(k, (str, type)):
                    raise TypeError(
                        "'kind' parameter must be a string/type or list of strings/types"
                    )
        elif not isinstance(kind, (str, type)):
            raise TypeError("'kind' parameter must be a string or type")

    if name is not None:
        if isinstance(name, (list, tuple)):
            for n in name:
                if not isinstance(n, str):
                    raise TypeError("'name' parameter must contain strings")
        elif not isinstance(name, str):
            raise TypeError(
                "'name' parameter must be a string or list/tuple of strings"
            )


def _matches_kind_pattern(element, kind_pattern):
    """Check if an element matches a kind pattern.

    Args:
        element: The element to check
        kind_pattern: String or type to match against

    Returns:
        bool: True if element matches the pattern
    """
    if isinstance(kind_pattern, str):
        # An element written as a Python subclass is still of its element kind, so match
        # the kind rather than the name of the user's class.
        if _matches_string(type(element).__name__, kind_pattern):
            return True
        return _matches_string(_element_kind(element), kind_pattern)
    elif isinstance(kind_pattern, type):
        return isinstance(element, kind_pattern)
    return False


def _element_kind(element):
    """The element kind, which a Python subclass keeps.

    Args:
        element: The element to inspect

    Returns:
        str: the name of the element type this is a kind of, e.g. ``"Drift"``
    """
    for base in type(element).__mro__:
        if getattr(elements, base.__name__, None) is base:
            return base.__name__
    return type(element).__name__


def _matches_name_pattern(element, name_pattern):
    """Check if an element matches a name pattern.

    Args:
        element: The element to check
        name_pattern: String pattern to match against

    Returns:
        bool: True if element matches the pattern
    """
    return (
        hasattr(element, "has_name")
        and element.has_name
        and _matches_string(element.name, name_pattern)
    )


def _check_element_match(element, kind, name):
    """Check if an element matches the given kind and name criteria.

    Args:
        element: The element to check
        kind: Kind criteria (str, type, list, tuple, or None)
        name: Name criteria (str, list, tuple, or None)

    Returns:
        bool: True if element matches any criteria (OR logic)
    """
    match = False

    # Check for 'kind' parameter
    if kind is not None:
        if isinstance(kind, (str, type)):
            # Single kind
            if _matches_kind_pattern(element, kind):
                match = True
        elif isinstance(kind, (list, tuple)):
            # Multiple kinds (OR logic)
            for k in kind:
                if _matches_kind_pattern(element, k):
                    match = True
                    break

    # Check for 'name' parameter (only if kind didn't match - OR logic)
    if name is not None and not match:
        if isinstance(name, str):
            # Single name
            if _matches_name_pattern(element, name):
                match = True
        elif isinstance(name, (list, tuple)):
            # Multiple names (OR logic)
            for name_item in name:
                if _matches_name_pattern(element, name_item):
                    match = True
                    break

    return match


def select(
    self,
    *,
    kind=None,
    name=None,
) -> FilteredElementsList:
    r"""Filter elements by type and name with OR-based logic.

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
                 Examples: "quad1", r"quad\d+", ["quad1", "quad2"], [r"quad\d+", "bend1"]
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

         lattice.select(name=r"quad\d+")  # Get elements matching pattern
         lattice.select(name=[r"quad\d+", "bend1"])  # Mix regex and strings

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

    # Handle keyword arguments for filtering
    if kind is not None or name is not None:
        # Validate parameters
        _validate_select_parameters(kind, name)

        matching_indices = []

        for i, element in enumerate(self):
            if _check_element_match(element, kind, name):
                matching_indices.append(i)

        return FilteredElementsList(self, matching_indices)

    # If no filtering criteria provided, return all elements
    all_indices = list(range(len(self)))
    return FilteredElementsList(self, all_indices)


def get_kinds(self) -> list[type]:
    """Get all unique element kinds in the list.

    Returns:
        list[type]: List of unique element types (sorted by name).
    """
    kinds = set()
    for element in self:
        kinds.add(type(element))
    return sorted(list(kinds), key=lambda t: t.__name__)


def count_by_kind(self, kind_pattern) -> int:
    """Count elements of a specific kind.

    Args:
        kind_pattern: The element kind to count. Can be:
            - String name (e.g., "Drift", "Quad") - supports exact match
            - Regex pattern (e.g., r".*Quad") - supports pattern matching
            - Element type (e.g., elements.Drift) - supports exact type match

    Returns:
        int: Number of elements of the specified kind.
    """
    count = 0
    for element in self:
        if _matches_kind_pattern(element, kind_pattern):
            count += 1
    return count


def has_kind(self, kind_pattern) -> bool:
    """Check if list contains elements of a specific kind.

    Args:
        kind_pattern: The element kind to check for. Can be:
            - String name (e.g., "Drift", "Quad") - supports exact match
            - Regex pattern (e.g., r".*Quad") - supports pattern matching
            - Element type (e.g., elements.Drift) - supports exact type match

    Returns:
        bool: True if at least one element of the specified kind exists.
    """
    for element in self:
        if _matches_kind_pattern(element, kind_pattern):
            return True
    return False


import math


def _filter_kwargs(d: dict) -> dict:
    """Filter a to_dict() result into valid constructor kwargs.

    Removes 'type' (not a constructor argument) and 'ds' when zero
    (thin elements don't accept ds).

    Args:
        d: Dictionary from element.to_dict()

    Returns:
        dict: Filtered dictionary suitable for element constructor
    """
    return {k: v for k, v in d.items() if k != "type" and (k != "ds" or v != 0.0)}


def _rad2deg(radians: float) -> float:
    """Convert radians to degrees."""
    return radians * 180.0 / math.pi


def _element_from_dict(d: dict):
    """Create an element from its to_dict() representation.

    Args:
        d: Dictionary from element.to_dict(), must include 'type' key

    Returns:
        Element instance of the appropriate type

    Raises:
        KeyError: If 'type' key is missing
        AttributeError: If element type is not found in elements module
    """
    element_type = d["type"]
    kwargs = _filter_kwargs(d)
    element_class = getattr(elements, element_type)
    return element_class(**kwargs)


def to_dicts(self) -> list[dict]:
    """Serialize the lattice to a list of dictionaries.

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
    # simple implementation:
    # return [el.to_dict() for el in self]

    # work-around for ExactSbend, PlaneXYRot, PRot, ThinDipole .to_dict() returning
    # radians not degrees; see _element_to_dict
    return [_element_to_dict(el) for el in self]


def _format_value(v):
    """Format a value for Python code generation.

    Converts AMReX SmallMatrix objects to lists.
    Returns a tuple of (formatted_value, matrix_dims) where matrix_dims
    is (rows, cols) from the SmallMatrix type, or None otherwise.
    """
    if hasattr(v, "to_numpy") and hasattr(v, "row_size"):
        rows = v.row_size
        cols = v.column_size
        arr = v.to_numpy()
        # Column vectors (Nx1) use flat list, matrices use nested list
        if cols == 1:
            return arr.flatten().tolist(), (rows, cols)
        else:
            return arr.T.tolist(), (rows, cols)
    return v, None


def to_py(self) -> str:
    """Generate Python code that recreates this lattice.

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
    # Check if we need amrex import for matrix types
    needs_amrex = False
    for d in self.to_dicts():
        kwargs = _filter_kwargs(d)
        for v in kwargs.values():
            if hasattr(v, "to_numpy") and hasattr(v, "row_size"):
                needs_amrex = True
                break
        if needs_amrex:
            break

    lines = ["from impactx import elements"]
    if needs_amrex:
        lines.append("import amrex.space3d as amr")

    lines.extend(
        [
            "",
            "",
            "def get_lattice():",
            '    """Return the lattice as a KnownElementsList."""',
            "    lattice = elements.KnownElementsList()",
            "    lattice.extend([",
        ]
    )

    # match the SmallMatrix element type to the build's precision
    real_t = "float" if Config.precision == "SINGLE" else "double"

    for d in self.to_dicts():
        element_type = d["type"]
        kwargs = _filter_kwargs(d)
        formatted_parts = []

        for k, v in kwargs.items():
            formatted_v, matrix_dims = _format_value(v)
            if matrix_dims:
                rows, cols = matrix_dims
                matrix_cls = f"amr.SmallMatrix_{rows}x{cols}_F_SI1_{real_t}"
                formatted_parts.append(f"{k}={matrix_cls}({repr(formatted_v)})")
            else:
                formatted_parts.append(f"{k}={repr(formatted_v)}")

        args_str = ", ".join(formatted_parts)
        lines.append(f"        elements.{element_type}({args_str}),")

    lines.append("    ])")
    lines.append("    return lattice")
    lines.append("")

    return "\n".join(lines)


def from_dicts(self, dicts: list[dict]):
    """Load and append elements from a list of dictionaries.

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
    self.extend([_element_from_dict(d) for d in dicts])


# ---------------------------------------------------------------------------
# Element-wise value comparison (== and isclose) for lattice containers.
#
# Duck-typed across container types: a KnownElementsList compares equal to a
# plain Python list of elements or to a FilteredElementsList as long as the
# element sequences match. This mirrors the typical select() workflow where
# users compare a filtered view against a hand-built reference list.
#
# No value-based __hash__ is defined: containers are mutable, and elements
# inside are mutable. The default identity-based hash inherited from ``object``
# is preserved (KnownElementsList from pybind11; FilteredElementsList from
# Python) so the existing weak-reference tracking of selections keeps working.
# ---------------------------------------------------------------------------


def _lattice_eq(self, other):
    """Element-wise equality with any iterable of elements.

    Duck-typed across container types: a ``KnownElementsList`` compares
    equal to a plain Python ``list`` of elements or to a
    ``FilteredElementsList`` as long as their lengths match and every
    pair of elements compares equal under ``==``. Returns
    ``NotImplemented`` for non-iterable operands so Python's
    reflected-equality fallback applies. Mutable containers are
    deliberately unhashable (``__hash__ = None``), matching the Python
    ``list`` convention.
    """
    if not hasattr(other, "__len__") or not hasattr(other, "__iter__"):
        raise ValueError(f"{other} does not implement __len__ or __iter__")
    if len(self) != len(other):
        return False
    for a, b in zip(self, other):
        if not (a == b):
            return False
    return True


def _lattice_isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
    """Tolerant element-wise comparison with any iterable of elements.

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
    if not hasattr(other, "__len__") or not hasattr(other, "__iter__"):
        raise ValueError(f"{other} does not implement __len__ or __iter__")
    if len(self) != len(other):
        return False
    for a, b in zip(self, other):
        if hasattr(a, "isclose"):
            if not a.isclose(
                b, rtol=rtol, atol=atol, ignore_attributes=ignore_attributes
            ):
                return False
        elif a != b:
            return False
    return True


# Apply to FilteredElementsList immediately. The pybind11 KnownElementsList is
# patched inside register_KnownElementsList_extension() below, since that is
# the sole entry point that receives the bound class.
#
# We deliberately do NOT touch __hash__: keeping the inherited identity-based
# hash means two value-equal containers may hash differently, but containers are
# mutable and are not intended as dict/set keys for value-based deduplication.
FilteredElementsList.__eq__ = _lattice_eq
FilteredElementsList.isclose = _lattice_isclose


def _lattice_init(self, elements=None):
    """Create a lattice, optionally filled with elements.

    The elements are shared, not copied, so the caller keeps handles to the very
    elements the lattice holds. Constructing goes through ``extend`` for exactly that
    reason: the C++ constructor cannot see the object being constructed, and so cannot
    record which Python objects own the elements.
    """
    _lattice_init_cxx(self)
    if elements is not None:
        if isinstance(elements, _elements_module.KnownElementsList) or isinstance(
            elements, (list, tuple)
        ):
            self.extend(elements)
        else:
            self.append(elements)


def register_KnownElementsList_extension(kel):
    """KnownElementsList helper methods"""
    from ..plot.Survey import plot_survey

    # Construction from an iterable of elements; see _lattice_init.
    global _lattice_init_cxx, _elements_module
    _lattice_init_cxx = kel.__init__
    from .. import impactx_pybind as _ix

    _elements_module = _ix.elements
    kel.__init__ = _lattice_init

    # register member functions for KnownElementsList
    kel.from_pals = from_pals
    kel.load_file = load_file
    kel.plot_survey = plot_survey

    # Serialization methods
    kel.to_dicts = to_dicts
    kel.to_py = to_py
    kel.from_dicts = from_dicts

    # Enhanced element selection methods
    kel.select = select
    kel.get_kinds = get_kinds
    kel.count_by_kind = count_by_kind
    kel.has_kind = has_kind

    # Element-wise == and isclose() (duck-typed across container types).
    # No custom __hash__: the inherited identity-based hash is kept so that
    # KnownElementsList and FilteredElementsList behave the same as Python
    # lists for the user-visible ``==``/``isclose`` API.
    kel.__eq__ = _lattice_eq
    kel.isclose = _lattice_isclose
