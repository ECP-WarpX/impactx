.. _examples-cyclotron:

Cyclotron
=========

This demonstrates a simple cyclotron as published by Ernest O. Lawrence and M. Stanley Livingston, *The Production of High Speed Light Ions Without the Use of High Voltages*, Phys. Rev. **40**, 19 (1932).
`DOI: 10.1103/PhysRev.40.19 <https://doi.org/10.1103/PhysRev.40.19>`__


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_cyclotron.py`` or
* ImpactX **executable** using an input file: ``impactx input_cyclotron.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_cyclotron.py
          :language: python3
          :caption: You can copy this file from ``examples/cyclotron/run_cyclotron.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_cyclotron.in
          :language: ini
          :caption: You can copy this file from ``examples/cyclotron/input_cyclotron.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_cyclotron.py``

   .. literalinclude:: analysis_cyclotron.py
      :language: python3
      :caption: You can copy this file from ``examples/cyclotron/analysis_cyclotron.py``.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``plot_cyclotron.py``

   .. literalinclude:: plot_cyclotron.py
      :language: python3
      :caption: You can copy this file from ``examples/cyclotron/plot_cyclotron.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/180287840-8561f6fd-278f-4856-abd8-04fbdb78c8ff.png
   :alt: focusing, defocusing and preserved emittane in our cyclotron cell benchmark.

   cyclotron transversal beam width and emittance evolution

.. figure:: https://user-images.githubusercontent.com/1353258/180287845-eb0210a7-2500-4aa9-844c-67fb094329d3.png
   :alt: focusing, defocusing and phase space rotation in our cyclotron cell benchmark.

   cyclotron transversal beam width and phase space evolution
