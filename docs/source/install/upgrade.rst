.. _install-upgrade:

Upgrade Guide
=============

This guide covers the changes that need attention when moving an existing script or
build to a newer ImpactX version: changed behavior, replaced APIs, and raised
requirements. Newest release first.

26.10
-----

``sim.lattice`` holds elements
""""""""""""""""""""""""""""""

The lattice now holds the element objects it is given, the way a Python list holds its
items. Previously it stored a copy of each element.

.. code-block:: python

   q = elements.Quad(ds=0.3, k=2.0)
   sim.lattice.append(q)

   assert sim.lattice[0] is q
   q.k = 3.0                       # applies to tracking

**What to check in your scripts:** adding the same variable more than once now places one
element at several positions, where it used to place independent copies. Retuning it
afterwards affects every one of those positions.

.. code-block:: python

   q = elements.Quad(ds=0.3, k=2.0)
   sim.lattice.extend([q, q])      # one element at two positions

The same holds for ``extend()``, which holds every element of the list it is given. A
repeated cell is the pattern most likely to be affected:

.. code-block:: python

   cell = [elements.Drift(ds=0.25), elements.Quad(ds=1.0, k=1.0)]
   sim.lattice.extend(cell * 3)    # 6 positions, 2 elements

This is harmless for elements that carry no per-position state, which is the common case:
a ``BeamMonitor`` at the head and tail of a lattice, or a drift repeated through a FODO
cell, behave as before. It matters where each position is expected to be tuned separately.

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

``sim.tracking_element`` yields the element being tracked, so writing to it from a hook
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

The same applies to the free functions ``elements.reverse(element)`` and
``push(pc, element)``, which now act on the element passed in.

Changing the lattice during tracking
""""""""""""""""""""""""""""""""""""

Adding or removing elements from a tracking hook raises, because tracking is walking the
sequence at that moment. Changing an element's own parameters from a hook is unaffected
and remains the intended way to retune a lattice mid-run.

New in this release
"""""""""""""""""""

- ``element.copy()`` on every element type, for a new element with the same configuration.
- The array-valued element parameters -- Fourier coefficients, multipole coefficients and
  polygon vertices -- can be set after construction, via their properties or the paired
  setter (``set_coefficients()``, ``set_vertices()``).
- ``sim.lattice`` supports the rest of the list API: slices, negative indices, ``insert()``,
  ``remove()``, ``index()``, ``count()``, ``in`` and ``del``.
