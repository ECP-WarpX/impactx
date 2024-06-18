.. _examples-fodo-programmable:

FODO Cell, Programmable Element
===============================

This implements the same FODO cell as the :ref:`stable FODO cell example <examples-fodo>`.
However, in the example here we define *additional user-defined, custom elements* (:py:class:`impactx.elements.Programmable`) from the :ref:`ImpactX Python APIs <usage-python>`.

More generally, ImpactX exposes all data structures through `pyAMReX for adding additional computation <https://pyamrex.readthedocs.io/en/latest/usage/compute.html>`__, enabling rapid prototyping of new elements on both CPU and GPU.


Run
---

This example can be run as a **Python** script: ``python3 fodo_programmable.py``.
For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix this line with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: run_fodo_programmable.py
   :language: python3
   :caption: You can copy this file from ``examples/fodo_programmable/run_fodo_programmable.py``.



Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo.py``

   .. literalinclude:: analysis_fodo.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_programmable/analysis_fodo.py``.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``plot_fodo.py``

   .. literalinclude:: plot_fodo.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_programmable/plot_fodo.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/180287840-8561f6fd-278f-4856-abd8-04fbdb78c8ff.png
   :alt: focusing, defocusing and preserved emittance in our FODO cell benchmark.

   FODO transversal beam width and emittance evolution

.. figure:: https://user-images.githubusercontent.com/1353258/180287845-eb0210a7-2500-4aa9-844c-67fb094329d3.png
   :alt: focusing, defocusing and phase space rotation in our FODO cell benchmark.

   FODO transversal beam width and phase space evolution
