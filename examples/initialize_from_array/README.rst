.. _examples-from-array:

Initialize a Beam from Arrays
=============================

This example demonstrates how a beam can be initialized in ImpactX from array-like structures.
This allows various applications of interest,
such as using a beam from a different simulation,
initializing a beam from file,
or creating a custom distribution.
In this example, a custom beam is specified at fixed t and added to ImpactX inside a
``with beam.at_fixed_t():`` block, which transforms it to the fixed-s coordinates of ImpactX
when the block ends (see :ref:`usage-howto-python-particle-data-fixed-t`).
The custom beam is a ring in x-y,
with radius r=2 mm,
radial width :math:`\sigma_r = 5\ \mathrm{\mu m}`;
Gaussian in :math:`p_x` and :math:`p_y` with momentum width :math:`\sigma_p=10`;
and chirped in z-pz with bunch length :math:`\sigma_z=1` mm,
mean energy about 10 GeV, 1% uncorrelated energy spread, and z-pz covariance of -0.18.

In specifying the beam at fixed t and transforming to fixed s,
it is assumed that the local and global coordinate frames align.
That is, the beam transverse directions x and y are the global x and y directions
and the beam z/t direction is the global z direction.
The run script uses a utility function to express the global particle coordinates relative to the reference particle.
The analysis script uses further utilities that reproduce the s-to-t transformation
done internally in ImpactX, given this assumption that the beam and global coordinate systems align.
These utility functions are provided in the following script:

.. dropdown:: Script ``transformation_utilities.py``

   .. literalinclude:: transformation_utilities.py
      :language: python3
      :caption: You can copy this file from ``examples/initialize_from_array/transformation_utilities.py``.


Run
---

This example can **only** be run with **Python**:

* **Python** script: ``python3 run_from_array.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. attention::

   In MPI-parallel simulations, ``beam.add_n_particles(...)`` is local to the MPI rank, spatial locality does not matter.
   Thus, you can add particles at any MPI rank, e.g., equally chuncked up for perfect load balancing.

   You do NOT want to add the same unique particle at multiple MPI ranks.

   When ImpactX needs to sort particles spatially, it will redistribute them over MPI ranks automatically during tracking.

.. literalinclude:: run_from_array.py
    :language: python3
    :caption: You can copy this file from ``examples/initialize_from_array/run_from_array.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analyze_from_array.py``

   .. literalinclude:: analyze_from_array.py
      :language: python3
      :caption: You can copy this file from ``examples/initialize_from_array/analyze_from_array.py``.

This uses the following utilities to read ImpactX output:

.. dropdown:: Script ``impactx_read_utilities.py``

   .. literalinclude:: impactx_read_utilities.py
      :language: python3
      :caption: You can copy this file from ``examples/initialize_from_array/impactx_read_utilities.py``.

Visualize
---------

We run the following script to visualize the ImpactX output and confirm the beam is properly initialized:

.. dropdown:: Script ``visualize_from_array.py``

   .. literalinclude:: visualize_from_array.py
      :language: python3
      :caption: You can copy this file from ``examples/initialize_from_array/visualize_from_array.py``.

The resulting phase space snapshots are shown in the following figure:

.. figure:: https://gist.githubusercontent.com/RTSandberg/613d16def3d025f9415d348a58bddba6/raw/f8d0dfe47f2a75063845b748df4cb3fb9f3b38bb/phase_space.png
   :alt: [fig:custom_beam] Phase space snapshots of custom beam from arrays, showing the ring in x-y, Gaussian in px-py, and linear correlation in z-pz.

   [fig:custom_beam] Phase space snapshots of custom beam from arrays, showing the ring in x-y, Gaussian in px-py, and linear correlation in z-pz.
