.. _running-cpp-parameters:

Input Parameters
================

.. note::
   The AMReX parser (see :ref:`running-cpp-parameters-parser`) is used for the right-hand-side of all input parameters that consist of one or more integers or floats, so expressions like ``<species_name>.density_max = "2.+1."`` and/or using user-defined constants are accepted.

.. _running-cpp-parameters-overall:

Overall simulation parameters
-----------------------------

* ``max_step`` (``integer``)
    The number of PIC cycles to perform.

* ``stop_time`` (``float``; in seconds)
    The maximum physical time of the simulation. Can be provided instead of ``max_step``. If both
    ``max_step`` and ``stop_time`` are provided, both criteria are used and the simulation stops
    when the first criterion is hit.

* ``amrex.abort_on_out_of_gpu_memory``  (``0`` or ``1``; default is ``1`` for true)
    When running on GPUs, memory that does not fit on the device will be automatically swapped to host memory when this option is set to ``0``.
    This will cause severe performance drops.
    Note that even with this set to ``1`` ImpactX will not catch all out-of-memory events yet when operating close to maximum device memory.
    `Please also see the documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters>`_.

.. _running-cpp-parameters-box:


Setting up the field mesh
-------------------------

* ``amr.n_cell`` (2 integers in 2D, 3 integers in 3D)
    The number of grid points along each direction (on the **coarsest level**)

* ``amr.max_level`` (``integer``, default: ``0``)
    When using mesh refinement, the number of refinement levels that will be used.

    Use 0 in order to disable mesh refinement.
    Note: currently, ``0`` and ``1`` are supported.

* ``amr.ref_ratio`` (`integer` per refined level, default: ``2``)
    When using mesh refinement, this is the refinement ratio per level.
    With this option, all directions are fined by the same ratio.

    Note: in development; currently, ``2`` is supported.

* ``amr.ref_ratio_vect`` (3 integers for x,y,z per refined level)
    When using mesh refinement, this can be used to set the refinement ratio per direction and level, relative to the previous level.

    Example: for three levels, a value of ``2 2 4 8 8 16`` refines the first level by 2-fold in x and y and 4-fold in z compared to the coarsest level (level 0/mother grid); compared to the first level, the second level is refined 8-fold in x and y and 16-fold in z.

    Note: in development; currently allowed value: ``2 2 2``.

* ``geometry.coord_sys`` (``integer``) optional (default ``0``)
    Coordinate system used by the simulation. 0 for Cartesian, 1 for cylindrical.

* ``geometry.prob_lo`` and ``geometry.prob_hi`` (2 floats in 2D, 3 floats in 3D; in meters)
    The extent of the full simulation box. This box is rectangular, and thus its
    extent is given here by the coordinates of the lower corner (``geometry.prob_lo``) and
    upper corner (``geometry.prob_hi``). The first axis of the coordinates is x
    (or r with cylindrical) and the last is z.


.. _running-cpp-parameters-bc:

Domain Boundary Conditions
--------------------------

.. note::

   TODO :-)

Lattice Elements
----------------

* ``lattice.elements`` (``list of strings``) optional (default: no elements)
    A list of names (one name per lattice element), in the order that they
    appear in the lattice.

* ``<element_name>.type`` (``string``)
    Indicates the element type for this lattice element. This should be one of:

        * ``drift`` for free drift. This requires this additional parameter:

            * ``<element_name>.ds`` (``float``, in meters) the segment length

        * ``quad`` for a quadrupole. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length

            * ``<element_name>.k`` (``float``, in inverse meters) the quadrupole strength

        * ``sbend`` for a bending magnet. This requires these additional parameters:

            * ``<element_name>.ds`` (``float``, in meters) the segment length

            * ``<element_name>.rc`` (``float``, in meters) the bend radius

.. _running-cpp-parameters-parallelization:

Distribution across MPI ranks and parallelization
-------------------------------------------------

* ``amr.max_grid_size`` (``integer``) optional (default ``128`)
    Maximum allowable size of each **subdomain**
    (expressed in number of grid points, in each direction).
    Each subdomain has its own ghost cells, and can be handled by a
    different MPI rank ; several OpenMP threads can work simultaneously on the
    same subdomain.

    If ``max_grid_size`` is such that the total number of subdomains is
    **larger** that the number of MPI ranks used, than some MPI ranks
    will handle several subdomains, thereby providing additional flexibility
    for **load balancing**.

    When using mesh refinement, this number applies to the subdomains
    of the coarsest level, but also to any of the finer level.


.. _running-cpp-parameters-parser:

Math parser and user-defined constants
--------------------------------------

ImpactX uses AMReX's math parser that reads expressions in the input file.
It can be used in all input parameters that consist of one or more integers or floats.
Integer input expecting boolean, 0 or 1, are not parsed.
Note that when multiple values are expected, the expressions are space delimited.
For integer input values, the expressions are evaluated as real numbers and the final result rounded to the nearest integer.
See `this section <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parser>`_ of the AMReX documentation for a complete list of functions supported by the math parser.

ImpactX constants
^^^^^^^^^^^^^^^^^

ImpactX will provide a few pre-defined constants, that can be used for any parameter that consists of one or more floats.

.. note::

   Develop, such as:

   ======== ===================
   q_e      elementary charge
   m_e      electron mass
   m_p      proton mass
   m_u      unified atomic mass unit (Dalton)
   epsilon0 vacuum permittivity
   mu0      vacuum permeability
   clight   speed of light
   pi       math constant pi
   ======== ===================

   See in WarpX the file ``Source/Utils/WarpXConst.H`` for the values.

User-defined constants
^^^^^^^^^^^^^^^^^^^^^^

Users can define their own constants in the input file.
These constants can be used for any parameter that consists of one or more integers or floats.
User-defined constant names can contain only letters, numbers and the character ``_``.
The name of each constant has to begin with a letter. The following names are used
by ImpactX, and cannot be used as user-defined constants: ``x``, ``y``, ``z``, ``X``, ``Y``, ``t``.
The values of the constants can include the predefined ImpactX constants listed above as well as other user-defined constants.
For example:

* ``my_constants.a0 = 3.0``
* ``my_constants.z_plateau = 150.e-6``
* ``my_constants.n0 = 1.e22``
* ``my_constants.wp = sqrt(n0*q_e**2/(epsilon0*m_e))``

Coordinates
^^^^^^^^^^^

Besides, for profiles that depend on spatial coordinates (the plasma momentum distribution or the laser field, see below ``Particle initialization`` and ``Laser initialization``), the parser will interpret some variables as spatial coordinates.
These are specified in the input parameter, i.e., ``density_function(x,y,z)`` and ``field_function(X,Y,t)``.

The parser reads python-style expressions between double quotes, for instance
``"a0*x**2 * (1-y*1.e2) * (x>0)"`` is a valid expression where ``a0`` is a
user-defined constant (see above) and ``x`` and ``y`` are spatial coordinates. The names are case sensitive. The factor
``(x>0)`` is ``1`` where ``x>0`` and ``0`` where ``x<=0``. It allows the user to
define functions by intervals.
Alternatively the expression above can be written as ``if(x>0, a0*x**2 * (1-y*1.e2), 0)``.

.. _running-cpp-parameters-particle:

Particle initialization
-----------------------

.. note::

   TODO :-)


.. _running-cpp-parameters-numerics:

Numerics and algorithms
-----------------------

.. note::

   TODO :-)


.. _running-cpp-parameters-diagnostics:

Diagnostics and output
----------------------

.. _running-cpp-parameters-diagnostics-insitu:

In-situ visualization
^^^^^^^^^^^^^^^^^^^^^

.. note::

   TODO :-)

.. _running-cpp-parameters-diagnostics-full:

.. note::

   TODO :-)

.. _running-cpp-parameters-diagnostics-reduced:

Reduced Diagnostics
^^^^^^^^^^^^^^^^^^^

``ReducedDiags`` allow the user to compute some reduced quantity (particle temperature, max of a field) and write a small amount of data to text files.

.. note::

   TODO :-)


.. _running-cpp-parameters-cp-restart:

Checkpoints and restart
-----------------------

.. note::

   ImpactX will support checkpoints/restart via AMReX.
   The checkpoint capability can be turned with regular diagnostics: ``<diag_name>.format = checkpoint``.

   * ``amr.restart`` (`string`)
       Name of the checkpoint file to restart from. Returns an error if the folder does not exist
       or if it is not properly formatted.

Intervals parser
----------------

.. note::

   TODO :-)

ImpactX can parse time step interval expressions of the form ``start:stop:period``, e.g.
``1:2:3, 4::, 5:6, :, ::10``.
A comma is used as a separator between groups of intervals, which we call slices.
The resulting time steps are the `union set <https://en.wikipedia.org/wiki/Union_(set_theory)>`_ of all given slices.
White spaces are ignored.
A single slice can have 0, 1 or 2 colons ``:``, just as `numpy slices <https://numpy.org/doc/stable/reference/generated/numpy.s_.html>`_, but with inclusive upper bound for ``stop``.

* For 0 colon the given value is the period

* For 1 colon the given string is of the type ``start:stop``

* For 2 colons the given string is of the type ``start:stop:period``

Any value that is not given is set to default.
Default is ``0`` for the start, ``std::numeric_limits<int>::max()`` for the stop and ``1`` for the
period.
For the 1 and 2 colon syntax, actually having values in the string is optional
(this means that ``::5``, ``100 ::10`` and ``100 :`` are all valid syntaxes).

All values can be expressions that will be parsed in the same way as other integer input parameters.

**Examples**

* ``something_intervals = 50`` -> do something at timesteps 0, 50, 100, 150, etc.
  (equivalent to ``something_intervals = ::50``)

* ``something_intervals = 300:600:100`` -> do something at timesteps 300, 400, 500 and 600.

* ``something_intervals = 300::50`` -> do something at timesteps 300, 350, 400, 450, etc.

* ``something_intervals = 105:108,205:208`` -> do something at timesteps 105, 106, 107, 108,
  205, 206, 207 and 208. (equivalent to ``something_intervals = 105 : 108 : , 205 : 208 :``)

* ``something_intervals = :`` or  ``something_intervals = ::`` -> do something at every timestep.

* ``something_intervals = 167:167,253:253,275:425:50`` do something at timesteps 167, 253, 275,
  325, 375 and 425.

This is essentially the python slicing syntax except that the stop is inclusive
(``0:100`` contains 100) and that no colon means that the given value is the period.

Note that if a given period is zero or negative, the corresponding slice is disregarded.
For example, ``something_intervals = -1`` deactivates ``something`` and
``something_intervals = ::-1,100:1000:25`` is equivalent to ``something_intervals = 100:1000:25``.
