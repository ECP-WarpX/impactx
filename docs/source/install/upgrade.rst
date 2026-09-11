.. _install-upgrade:

Upgrade Guide
=============

This guide covers the changes that need attention when moving an existing script or
build to a newer ImpactX version: changed behavior, replaced APIs, and raised
requirements. Newest release first.

26.10
-----

:py:attr:`sim.lattice <impactx.ImpactX.lattice>` holds elements
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The lattice now holds the element objects it is given, the way a Python list holds its
items. Previously it stored a copy of each element.

.. code-block:: python

   q = elements.Quad(ds=0.3, k=2.0)
   sim.lattice.append(q)

   assert sim.lattice[0] is q
   q.k = 3.0                       # applies to tracking

**What to check in your scripts:** previously, when storing an element in a variable and adding that variable more than once used to place independent copies.
Now, this *places the same element at several positions*.
Changing the variable (tuning the element) afterwards affects every one of those positions.

.. code-block:: python

   q = elements.Quad(ds=0.3, k=2.0)
   sim.lattice.extend([q, q])      # one element at two positions

The same holds for :py:meth:`~impactx.elements.KnownElementsList.extend`, which holds
every element of the list it is given. A
repeated cell is the pattern most likely to be affected:

.. code-block:: python

   cell = [elements.Drift(ds=0.25), elements.Quad(ds=1.0, k=1.0)]
   sim.lattice.extend(cell * 3)    # 6 positions, 2 elements

This is harmless for elements that carry no per-position state, which is the common case:
a :py:class:`~impactx.elements.BeamMonitor` at the head and tail of a lattice, or a drift
repeated through a FODO cell, behave as before. It matters where each position is expected to be tuned separately.

For separate elements, construct one per position or add a copy:

.. code-block:: python

   sim.lattice.append(elements.Quad(ds=0.3, k=2.0))
   sim.lattice.append(elements.Quad(ds=0.3, k=2.0))

   template = elements.Quad(ds=0.3, k=2.0)
   sim.lattice.extend([template.copy(), template.copy()])

To repeat a cell as independent elements, copy the list per repetition:

.. code-block:: python

   for _ in range(3):
       sim.lattice.extend([element.copy() for element in cell])

``sim.tracking_element`` writes take effect
"""""""""""""""""""""""""""""""""""""""""""

:py:attr:`sim.tracking_element <impactx.ImpactX.tracking_element>` yields the element being tracked, so writing to it from a hook
changes the push that follows. It previously yielded a copy, and such writes were
discarded.

If a hook retunes an element relative to its own current value, and that element occupies
several positions, read the reference value once before tracking instead of reading it
back from the element:

.. code-block:: python

   phase_shift = rf.phase          # once, before tracking

   def hook_before_element(sim):
       element = sim.tracking_element
       if type(element) is elements.RFCavity:
           element.phase = optimize(sim.beam.ref, element) + phase_shift

The same applies to the free functions :py:func:`impactx.reverse` and
:py:func:`impactx.push`, which now act on the element passed in.

Changing the lattice during tracking
""""""""""""""""""""""""""""""""""""

Adding or removing elements from a tracking hook raises, because tracking is walking the
sequence at that moment. Changing an element's own parameters from a hook is unaffected
and remains the intended way to retune a lattice mid-run.

``sim.lattice`` needs its simulation
""""""""""""""""""""""""""""""""""""

A lattice reached through :py:attr:`sim.lattice <impactx.ImpactX.lattice>` belongs to that simulation and does not keep it
alive. Using the lattice after the simulation is gone raises, where it used to keep
working:

.. code-block:: python

   def build():
       sim = ImpactX()
       sim.lattice.extend(cell)
       return sim.lattice        # the simulation ends here

   lattice = build()
   len(lattice)                  # raises

Return the simulation, or keep a reference to it for as long as you use its lattice.
Alternatively, keep the lattice alive across simulations by storing it in a pure Python
list and call :py:meth:`~impactx.elements.KnownElementsList.extend` directly before tracking, see
:ref:`usage-howto-lattice-manipulation`.

Selections describe the lattice they were taken from
""""""""""""""""""""""""""""""""""""""""""""""""""""

:py:meth:`~impactx.elements.KnownElementsList.select` returns a set of positions. Adding, removing or moving elements afterwards
makes those positions describe something else, so the selection raises rather than
reporting the wrong elements:

.. code-block:: python

   quads = sim.lattice.select(kind="Quad")
   sim.lattice.append(elements.Drift(ds=0.1))
   quads[0]                      # raises; select() again

Changing an element's own parameters moves nothing and leaves a selection usable.

``select(kind=...)`` matches the element kind
"""""""""""""""""""""""""""""""""""""""""""""

An element written as a Python subclass is still of its element kind, so
:py:meth:`~impactx.elements.KnownElementsList.select` with ``kind="Drift"`` or ``kind=elements.Drift`` now matches it. The name of the
subclass matches as well.

``replace_each()`` needs a template it can copy
"""""""""""""""""""""""""""""""""""""""""""""""

:py:meth:`~impactx.elements.FilteredElementsList.replace_each` copies the template into
each selected position, so a Python subclass must define
:py:meth:`~impactx.elements.Element.copy`. Passing one that does not raises, where it previously inserted plain
base-class elements.

``finalize()`` always empties the lattice
"""""""""""""""""""""""""""""""""""""""""

``sim.finalize()`` releases the elements whether or not
:py:meth:`~impactx.ImpactX.init_grids` was ever called.
They have been finalized by then, so keeping them would leave elements that are done with,
e.g., a :py:class:`~impactx.elements.BeamMonitor` with its output closed, still in the
lattice.

New in this release
"""""""""""""""""""

- :py:meth:`~impactx.elements.Element.copy` on every element type, for a new element with
  the same configuration.
  Keyword arguments give the copy a different value for a parameter, so that one element
  can serve as the template for many: ``[quad.copy(k=k) for k in scan]``.
- The array-valued element parameters (Fourier coefficients, multipole coefficients and
  polygon vertices) can be set after construction, via their properties or the paired
  setter (``set_coefficients()``, ``set_vertices()``).
- ``sim.lattice`` supports the rest of the list API: slices, negative indices,
  :py:meth:`~impactx.elements.KnownElementsList.insert`, :py:meth:`~impactx.elements.KnownElementsList.remove`,
  :py:meth:`~impactx.elements.KnownElementsList.index`, :py:meth:`~impactx.elements.KnownElementsList.count`,
  :py:meth:`in <impactx.elements.KnownElementsList.__contains__>`,
  :py:meth:`del <impactx.elements.KnownElementsList.__delitem__>`,
  :py:meth:`lattice[i] = element <impactx.elements.KnownElementsList.__setitem__>` and
  :py:meth:`reversed() <impactx.elements.KnownElementsList.__reversed__>`.
  :py:meth:`~impactx.elements.KnownElementsList.pop_back` returns the element it removed, and a position that is not
  there raises ``IndexError``.
- ``lattice.generation`` counts the structural edits made to a lattice.
