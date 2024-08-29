.. _examples-chicane:

Chicane
=======

This is the Berlin-Zeuthen magnetic bunch compression chicane, which is a standardized community benchmark.

`All parameters can be found online <https://www.desy.de/csr/>`__.
A 5 GeV electron bunch with normalized transverse rms emittance of 1 um undergoes longitudinal compression by a factor of 10 in a standard 4-bend chicane.

The emittances should be preserved, and the rms pulse length should decrease by the compression factor (10).

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

We also have :ref:`a variation of this test that includes CSR effects in the bending magnets <examples-chicane-csr>`.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_chicane.py`` or
* ImpactX **executable** using an input file: ``impactx input_chicane.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_chicane.py
          :language: python3
          :caption: You can copy this file from ``examples/chicane/run_chicane.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_chicane.in
          :language: ini
          :caption: You can copy this file from ``examples/chicane/input_chicane.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_chicane.py``

   .. literalinclude:: analysis_chicane.py
      :language: python3
      :caption: You can copy this file from ``examples/chicane/analysis_chicane.py``.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``plot_chicane.py``

   .. literalinclude:: plot_chicane.py
      :language: python3
      :caption: You can copy this file from ``examples/chicane/plot_chicane.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/180332191-f9ce11fc-8c56-4713-a91a-2ad12ab09805.png
   :alt: Chicane floorplan, beam width and restored emittane in our Chicane benchmark

   (top) Chicane floorplan.
   (bottom) Chicane beam width and emittance evolution.

.. figure:: https://user-images.githubusercontent.com/1353258/181611473-754dde72-3281-453b-9d9a-43317a5a49f2.png
   :alt: Beam transversal compression in our chicane example.

   Chicane beam width and emittance evolution
