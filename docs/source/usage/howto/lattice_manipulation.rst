.. _usage-howto-lattice-manipulation:

Manipulate a Lattice
====================

:py:attr:`sim.lattice <impactx.ImpactX.lattice>` is a sequence of elements that behaves like a Python list.
When you add elements to ``sim.lattice`` it will **hold references to** the Python variables it is given rather than copies of them.
Concretely: if you create an element as a Python variable, append it to ``sim.lattice``, you can continue to change and manipulate the element that now sits inside ``sim.lattice`` through the variable (see below).
_The same element (variable) may sit at several lattice positions_, making tuning and ramping operations for whole channels, ring arcs, etc. easy.

This page collects the operations that come up most often:
editing a lattice in place, building one lattice per run, and keeping a lattice for longer
than the simulation that used it.
For the full signature of every method, see
:py:class:`impactx.elements.KnownElementsList`.

Everything is a Reference
-------------------------

Adding an element borrows the object to the lattice.
Changing it afterwards changes what is tracked:

.. code-block:: python

   q = elements.Quad(ds=0.3, k=2.0)
   sim.lattice.append(q)

   q.k = 4.0              # sim.lattice references variable q, this changes the lattice

   sim.track_particles()  # tracks through Quad(ds=0.3, k=4.0)

Adding the same variable twice therefore places **one** element at two positions, and
retuning it changes both:

.. code-block:: python

   sim.lattice.extend([q, q])    # one element, two positions

For elements that are tuned separately, construct one per position or add a copy:

.. code-block:: python

   sim.lattice.append(elements.Quad(ds=0.3, k=2.0))
   sim.lattice.append(elements.Quad(ds=0.3, k=2.0))   # two distinct elements

   template = elements.Quad(ds=0.3, k=2.0)
   sim.lattice.extend([template.copy(), template.copy()])  # also two distinct elements

Repeating a Cell
----------------

:py:meth:`~impactx.elements.KnownElementsList.extend` holds every element of the list it is given, so multiplying a list
repeats its elements rather than duplicating them:

.. code-block:: python

   cell = [elements.Drift(ds=0.25), elements.Quad(ds=1.0, k=1.0)]

   sim.lattice.extend(cell * 3)          # 6 positions, 2 elements used repeatedly
   sim.lattice[0].ds = 0.5               # changes all three cells

That is usually what a periodic channel wants: retuning the cell retunes the channel.
To tune the cells separately, copy the elements once per repetition:

.. code-block:: python

   for _ in range(3):
       sim.lattice.extend([element.copy() for element in cell])

.. note::

   Copying the *list* does not copy the elements in it.
   ``cell.copy()``, ``cell[:]`` and ``list(cell)`` each give a new list holding the same
   element objects, so this still places the very elements ``cell`` holds:

   .. code-block:: python

      sim.lattice.extend(cell.copy())   # a new list, the same elements

   `copy.copy() <https://docs.python.org/3/library/copy.html#copy.copy>`__ also copies
   the list, but references the same elements.
   Copy the elements, as above, to get independent ones.

Editing a Lattice in Place
--------------------------

The list operations work as they do on a ``list``, and act on the lattice being tracked:

.. code-block:: python

   sim.lattice[1].k = 2.5                                   # retune one element
   sim.lattice.insert(0, elements.Marker(name="start"))     # add at a position
   del sim.lattice[-1]                                      # remove the last
   sim.lattice[2:4] = [elements.Drift(ds=0.75)]             # replace a range

   sim.lattice[3] = elements.Quad(ds=1.0, k=-1.0)           # swap one element out

Indexing, slicing, iteration, ``len()``, :py:meth:`in <impactx.elements.KnownElementsList.__contains__>`,
:py:meth:`~impactx.elements.KnownElementsList.insert`, :py:meth:`~impactx.elements.KnownElementsList.remove`, :py:meth:`~impactx.elements.KnownElementsList.index`,
:py:meth:`~impactx.elements.KnownElementsList.count`, :py:meth:`del <impactx.elements.KnownElementsList.__delitem__>` and
:py:meth:`reversed() <impactx.elements.KnownElementsList.__reversed__>` are all available.
Membership, :py:meth:`~impactx.elements.KnownElementsList.index`, :py:meth:`~impactx.elements.KnownElementsList.count` and :py:meth:`~impactx.elements.KnownElementsList.remove` match on
the element itself, so an element that merely has the same parameters as another is not
mistaken for it.

To work on a group of elements at once, select them with
:py:meth:`~impactx.elements.KnownElementsList.select`, which gives a
:py:class:`impactx.elements.FilteredElementsList`:

.. code-block:: python

   sim.lattice.select(kind="Quad").replace_each(elements.Drift(ds=1.0))
   sim.lattice.select(name="corrector.*").delete()

A name is matched exactly, or as a regular expression when it looks like one, so
``"corrector.*"`` selects every corrector and ``"corrector1"`` selects that one.

A selection is a set of positions, so it describes the lattice as it was when it was taken.
Adding, removing or moving elements afterwards makes it stale, and using it then raises;
take a new one with :py:meth:`~impactx.elements.KnownElementsList.select`.
Changing an element's own parameters moves nothing and leaves a selection usable.

.. _usage-howto-lattice-manipulation-scan:

One Lattice per Run
-------------------

A parameter scan wants each run to differ in one knob and share everything else.
The simplest way is a function that builds a fresh cell per run:

.. code-block:: python

   def fodo_cell(k):
       return [
           elements.Drift(ds=0.25, name="d1"),
           elements.Quad(ds=1.0, k=k, name="q1"),
           elements.Drift(ds=0.5, name="d2"),
           elements.Quad(ds=1.0, k=-k, name="q2"),
           elements.Drift(ds=0.25, name="d3"),
       ]

   for k in [0.8, 0.9, 1.0, 1.1]:
       sim = ImpactX()
       # ... set up the beam ...
       sim.lattice.extend(fodo_cell(k))
       sim.track_particles()
       sim.finalize()

Every run gets its own elements, so nothing a run does can reach the runs around it.

When the elements already exist -- loaded from a lattice file, or built once at the top of
the script -- derive each run from those instead.
Copy the elements the scan changes and share the rest:

.. code-block:: python

   d1 = elements.Drift(ds=0.25, name="d1")
   q1 = elements.Quad(ds=1.0, k=1.0, name="q1")
   d2 = elements.Drift(ds=0.5, name="d2")
   q2 = elements.Quad(ds=1.0, k=-1.0, name="q2")
   d3 = elements.Drift(ds=0.25, name="d3")

   for k in [0.8, 0.9, 1.0, 1.1]:
       sim = ImpactX()
       # ... set up the beam ...
       sim.lattice.extend([d1, q1.copy(k=k), d2, q2.copy(k=-k), d3])
       sim.track_particles()
       sim.finalize()

The quadrupoles are copied because each run gives them a different strength; the drifts are
not, because no run changes them.
``q1`` and ``q2`` keep their original values throughout, so the scan reads off the loop
variable rather than whatever the previous run left behind.

Keyword arguments to :py:meth:`~impactx.elements.Element.copy` say how the copy differs
from the element it was made from, so that one element supplies the values for many:

.. code-block:: python

   scan = [q1.copy(k=k) for k in [0.8, 0.9, 1.0, 1.1]]

   rf1 = elements.RFCavity(name="rf1", ds=1.3, escale=20.0, freq=1.3e9, phase=0.0,
                           cos_coefficients=[2.0], sin_coefficients=[0.0])
   cavities = [rf1] + [rf1.copy(name=f"rf{i}") for i in range(2, 5)]

Parameters that only mean something as a pair are given together, as the constructor takes
them:

.. code-block:: python

   square = elements.PolygonAperture(vertices_x=[-1.0, 1.0, 1.0, -1.0, -1.0],
                                     vertices_y=[-1.0, -1.0, 1.0, 1.0, -1.0])
   wider = square.copy(vertices_x=[-2.0, 2.0, 2.0, -2.0, -2.0],
                       vertices_y=[-2.0, -2.0, 2.0, 2.0, -2.0])

Setting a parameter the element does not have raises an error, e.g., for typos.

.. note::

   :py:meth:`element.copy() <impactx.elements.Element.copy>` gives a new element with the same configuration.
   A copy of a beam monitor does not inherit an already-open output file, and an element
   written as a Python subclass must define ``copy()`` for this to work.

Keeping a Lattice Longer Than a Simulation
------------------------------------------

A lattice accessed through :py:attr:`sim.lattice <impactx.ImpactX.lattice>` belongs to that
simulation, using it after the simulation is gone raises an error.
To reuse a lattice across simulations, keep the elements in a plain Python list and extend
each simulation's lattice from it before tracking:

.. code-block:: python

   cell = [                      # a plain Python list, owned by your script
       elements.Drift(ds=0.25, name="d1"),
       elements.Quad(ds=1.0, k=1.0, name="q1"),
       elements.Drift(ds=0.25, name="d2"),
   ]

   for run in range(runs):
       sim = ImpactX()
       # ... set up the beam ...
       sim.lattice.extend(cell)  # the same elements each run
       sim.track_particles()
       sim.finalize()

The list outlives every simulation, and each run tracks the same element objects.
Retuning a lattice list between runs therefore applies to the runs that follow.
Use :py:meth:`element.copy() <impactx.elements.Element.copy>` where the runs
must be independent, as in the scan above.

.. note::

   Returning :py:attr:`sim.lattice <impactx.ImpactX.lattice>` from a function does not work,
   because the simulation ends
   with the function:

   .. code-block:: python

      def build():
          sim = ImpactX()
          sim.lattice.extend(cell)
          return sim.lattice     # the simulation ends here

      lattice = build()
      len(lattice)               # raises

   Return the simulation as well, or keep the elements in your own list as above.
