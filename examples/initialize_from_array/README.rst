.. _examples-from-array:

Initialize a beam from arrays
=============================

This example demonstrates how a beam can be initalized in ImpactX from array-like structures.
This allows various applications of interest,
such as using a beam from a different simulation,
initializing a beam from file,
or creating a custom distribution.
This example includes a set of utilities for transforming the beam to the fixed-s coordinates of ImpactX.

In this example, a custom beam is specified at fixed t, transformed to fixed s, and
then loaded in ImpactX. 
The custom beam is a ring in x-y,
with radius r=2 mm,
radial width :math:`\sigma_r = 5\ \mu`m;
Gaussian in :math:`p_x` and :math:`p_y` with momentum width :math:`\sigma_p=10`;
and chirped in z-pz with bunch length :math:`\sigma_z=1` mm,
mean energy about 10 GeV, 1% uncorrelated energy spread, and z-pz covariance of -0.18.

In specifying the beam at fixed t and transforming to fixed s,
it is assumed that the local and global coordinate frames align.
That is, the beam transverse directions x and y are the global x and y directions
and the beam z/t direction is the global z direction.
The transformation utility function reproduces the t-to-s and s-to-t transformations
done internally in ImpactX given this assumption that the beam and global coordinate systems align.
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
