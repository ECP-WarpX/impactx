.. _usage-workflows-add-element:

Add New Beamline Elements
=========================

In ImpactX, one can easily add new beamline elements as a user.
There are multiple ways to add new elements to ImpactX, you can pick the one that fits your needs best.

The workflows described here apply both for thin kicks or thick elements.
Thick elements can also use soft-edged fringe fields (see `existing soft-edged elements for implementation details <https://github.com/ECP-WarpX/impactx/tree/development/src/particles/elements>`__).


Python Programmable Element
---------------------------

Using the :ref:`ImpactX Python interface <usage-python>`, a custom element named :py:class:`impactx.elements.Programmable` can be defined to advance particles using NumPy, CuPy, Numba, PyTorch or any other compatible Python library.

The Programmable element can implement a custom element in two ways:

* Push the whole container, by assigning a ``push`` function or
* Push the reference particle and beam particles in two individual functions (``beam_particles`` and ``ref_particle``).

Per ImpactX convention, the reference particle is updated *before* the beam particles are pushed.

Detailed examples that show usage of the programmable element are:

* :ref:`FODO cell <examples-fodo-programmable>`: implements a user-defined drift
* :ref:`15 stage laser-plasma accelerator <examples-ml-surrogate>`: implements a user-defined LPA accelerator element using a neural network surrogate via PyTorch

Detailed particle computing interfaces are presented in the `pyAMReX examples <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__.


Linear Map
----------

.. note::

   We plan to add a simple, linear map element that can be configured in user input.
   Follow `issue #538 <https://github.com/ECP-WarpX/impactx/issues/538>`__ for progress.


C++ Element
-----------

Adding a new beamline element directly to the C++ code base of ImpactX is straight forward and described in the following.

We store all beamline elements under `src/particles/elements/ <https://github.com/ECP-WarpX/impactx/tree/development/src/particles/elements>`__.

Let's take a look at an example, the `Drift <https://impactx.readthedocs.io/en/latest/_static/doxyhtml/structimpactx_1_1_drift.html>`__ implementation.
To simplify the logic, we use so-called `mixin classes <https://en.wikipedia.org/wiki/Mixin>`__, which provide commonly used logic for `parallelization, thin/thick elements, alignment error support, etc <https://impactx.readthedocs.io/en/latest/_static/doxyhtml/namespaceimpactx_1_1elements.html>`__.

.. literalinclude:: ../../../../src/particles/elements/Drift.H
   :language: cpp
   :dedent: 4
   :start-at: struct Drift
   :end-at: static constexpr auto type = "Drift";

After this brief boilerplate, our beamline elements implement three simple parts:

* a constructor: storing element options
* a single-particle operator: pushing the beam particles
* a reference-particle operator: pushing the reference particle

.. dropdown:: Example Element: Drift.H
   :color: light
   :icon: info
   :animate: fade-in-slide-down

   .. literalinclude:: ../../../../src/particles/elements/Drift.H
      :language: cpp

As a last step, we expose our C++ beamline elements to Python in `src/python/elements.cpp <https://github.com/ECP-WarpX/impactx/blob/development/src/python/elements.cpp>`__.

.. dropdown:: Python Binding: Drift
   :color: light
   :icon: info
   :animate: fade-in-slide-down

   .. literalinclude:: ../../../../src/python/elements.cpp
      :language: cpp
      :dedent: 8
      :start-at: py::class_<Drift, elements::Thick, elements::Alignment> py_Drift(me, "Drift");
      :end-at: register_beamoptics_push(py_Drift);

Pull requests that added a new element and can be taken as examples are:

  * `Chromatic Plasma Lens <https://github.com/ECP-WarpX/impactx/pull/514>`__
  * `Thin-Kick Dipole <https://github.com/ECP-WarpX/impactx/pull/472>`__
  * `Chromatic Elements for Drift, Quad, Uniform Focusing+Solenoid <https://github.com/ECP-WarpX/impactx/pull/356>`__
  * `Quadrupole with Soft-Edge Fringe Fields <https://github.com/ECP-WarpX/impactx/pull/322>`__
  * other pull requests under the `component: elements <https://github.com/ECP-WarpX/impactx/pulls?q=is%3Apr+label%3A%22component%3A+elements%22+is%3Amerged+>`__ label
