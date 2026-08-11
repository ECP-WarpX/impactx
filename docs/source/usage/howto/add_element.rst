.. _usage-howto-add-element:

Add New Beamline Elements
=========================

In ImpactX, one can easily add new beamline elements as a user.
There are multiple ways to add new elements to ImpactX, you can pick the one that fits your needs best.

The workflows described here apply both for thin kicks or thick elements.
Thick elements can also use soft-edged fringe fields (see `existing soft-edged elements for implementation details <https://github.com/BLAST-ImpactX/impactx/tree/development/src/elements>`__).


.. _usage-howto-add-element-linmap:

Linear Map
----------

A custom linear element can be provided by specifying the 6x6 linear transport matrix :math:`R` as an input.
See the :ref:`FODO cell example <examples-fodo-userdef>` for a demonstration of the custom linear element.

The matrix elements :math:`R(i,j)` are indexed beginning with 1, so that :math:`i,j=1,2,3,4,5,6`.

The matrix :math:`R` multiplies the phase space vector :math:`(x,p_x,y,p_y,t,p_t)`, where coordinates :math:`(x,y,t)` have units of m and momenta :math:`(p_x,p_y,p_t)` are dimensionless.
So, for example, :math:`R(1,1)` is dimensionless, and :math:`R(1,2)` has units of m.


.. note::

   If a user-provided linear map is used, it is up to the user to ensure that the 6x6 transport matrix is symplectic.
   If a more general form of user-defined transport is needed, the :ref:`Python Programmable Element <usage-howto-add-element-python>` and the :ref:`C++ Element <usage-howto-add-element-cxx>` provide a more general approach.


.. _usage-howto-add-element-python:

Python Programmable Element
---------------------------

Using the :ref:`ImpactX Python interface <usage-python>`, a custom element named :py:class:`impactx.elements.Programmable` can be defined to advance particles using NumPy, CuPy, Numba, PyTorch or any other compatible Python library.

The Programmable element can implement a custom element in two ways:

* Push the whole container, by assigning a ``push`` function or
* Push the reference particle and beam particles in two individual functions (``beam_particles`` and ``ref_particle``).

Per ImpactX convention, the reference particle is updated *before* the beam particles are pushed.

.. note::

   The Programmable element is for *replacing* a beamline element's particle push.
   For in-situ **analysis** of the beam, use :py:attr:`sim.hook <impactx.ImpactX.hook>`
   callbacks with ``sim.beam`` instead — see :ref:`usage-howto-python-extend`.

Detailed examples that show usage of the programmable element are:

* :ref:`FODO cell <examples-fodo-programmable>`: implements a user-defined drift
* :ref:`15 stage laser-plasma accelerator <examples-ml-surrogate>`: implements a user-defined LPA accelerator element using a neural network surrogate via PyTorch

Detailed particle computing interfaces are presented in the `pyAMReX examples <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__.


.. _usage-howto-add-element-cxx:

C++ Element
-----------

Adding a new beamline element directly to the C++ code base of ImpactX is straight forward and described in the following.

We store all beamline elements under `src/elements/ <https://github.com/BLAST-ImpactX/impactx/tree/development/src/elements>`__.

Let's take a look at an example, the `Drift <https://impactx.readthedocs.io/en/latest/_static/doxyhtml/structimpactx_1_1elements_1_1_drift.html>`__ implementation.
To simplify the logic, we use so-called `mixin classes <https://en.wikipedia.org/wiki/Mixin>`__, which provide commonly used logic for `parallelization, thin/thick elements, alignment error support, etc <https://impactx.readthedocs.io/en/latest/_static/doxyhtml/namespaceimpactx_1_1elements_1_1mixin.html>`__.

.. literalinclude:: ../../../../src/elements/Drift.H
   :language: cpp
   :dedent: 4
   :start-at: struct Drift
   :end-at: static constexpr auto type = "Drift";

After this brief boilerplate, our beamline elements implement these parts:

#. a constructor: storing element options
#. for particle tracking:

   * a single-particle operator: pushing the beam particles
   * a reference-particle operator: pushing the reference particle
#. for envelope tracking: a linear transport map

.. dropdown:: Example Element: Drift.H
   :color: light
   :icon: info
   :animate: fade-in-slide-down

   .. literalinclude:: ../../../../src/elements/Drift.H
      :language: cpp

As a last step, we expose our C++ beamline elements to Python in `src/python/elements.cpp <https://github.com/BLAST-ImpactX/impactx/blob/development/src/python/elements.cpp>`__.

.. dropdown:: Python Binding: Drift
   :color: light
   :icon: info
   :animate: fade-in-slide-down

   .. literalinclude:: ../../../../src/python/elements.cpp
      :language: cpp
      :dedent: 4
      :start-at: py::class_<Drift, elements::mixin::Named, elements::mixin::Thick, elements::mixin::Alignment, elements::mixin::PipeAperture, py::smart_holder> py_Drift(me, "Drift");
      :end-at: register_push(py_Drift);

Pull requests that added a new element and can be taken as examples are:

* `Chromatic Plasma Lens <https://github.com/BLAST-ImpactX/impactx/pull/514>`__
* `Thin-Kick Dipole <https://github.com/BLAST-ImpactX/impactx/pull/472>`__
* `Chromatic Elements for Drift, Quad, Uniform Focusing+Solenoid <https://github.com/BLAST-ImpactX/impactx/pull/356>`__
* `Quadrupole with Soft-Edge Fringe Fields <https://github.com/BLAST-ImpactX/impactx/pull/322>`__
* other pull requests under the `component: elements <https://github.com/BLAST-ImpactX/impactx/pulls?q=is%3Apr+label%3A%22component%3A+elements%22+is%3Amerged+>`__ label
