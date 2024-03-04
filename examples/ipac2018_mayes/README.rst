.. _examples-gaussian-SCfields:

Space Charge Fields in a 3D Gaussian Bunch
===========================================

This example is based on Fig. 1 of:

C. E. Mayes, R. D. Ryne, and D. C. Sagan, "3D Space Charge in BMAD", in Proc. IPAC2018, Vancouver, BC, Canada,
doi:10.18429/JACoW-IPAC2018-THPAK085

This is a test of the ImpactX Poisson solver.  The space charge fields are computed within a 1 nC electron bunch with a 3D Gaussian charge distribution
(at rest) and compared with analytical results for several values of transverse:longitudinal beam aspect ratio (Fig. 1).

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run only as a:

* **Python** script: ``python3 run_fieldplots_ipac2018.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fieldplots_ipac2018.py
          :language: python3
          :caption: You can copy this file from ``examples/ipac2018_mayes/run_fieldplots_ipac2018.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fieldplots_ipac2018.py``

   .. literalinclude:: analysis_fieldplots_ipac2018.py
      :language: python3
      :caption: You can copy this file from ``examples/ipac2018_mayes/analysis_fieldplots_ipac2018.py``.



.. _examples-guassian-SCdrift:

Space Charge in a Coasting 3D Gaussian Bunch        
============================================= 

This example is based on Fig. 4 of:

C. E. Mayes, R. D. Ryne, and D. C. Sagan, "3D Space Charge in BMAD", in Proc. IPAC2018, Vancouver, BC, Canada,
doi:10.18429/JACoW-IPAC2018-THPAK085

A cold (zero momentum spread) 1 nC electron bunch with a 3D Gaussian distribution and a (total) reference energy of 10 MeV is allowed to expand in a drift of length 1 m.
The final phase space is compared against the expected result (Fig. 4).

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsil>


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_expanding_ipac2018.py`` or
* ImpactX **executable** using an input file: ``impactx input_expanding_ipac2018.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_expanding_ipac2018.py
          :language: python3
          :caption: You can copy this file from ``examples/ipac2018_mayes/run_expanding_ipac2018.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_expanding_ipac2018.in
          :language: ini
          :caption: You can copy this file from ``examples/ipac2018_mayes/input_expanding_ipac2018.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_expanding_ipac2018.py``

   .. literalinclude:: analysis_expanding_ipac2018.py
      :language: python3
      :caption: You can copy this file from ``examples/ipac2018_mayes/analysis_expanding_ipac2018.py``.


Visualize
---------

You can run the following script to visualize the final beam distribution:

.. dropdown:: Script ``plot_expanding_ipac2018.py``

   .. literalinclude:: plot_expanding_ipac2018.py
      :language: python3
      :caption: You can copy this file from ``examples/ipac2018_mayes/plot_expanding_ipac2018.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/294003440-b16185c7-2573-48d9-8998-17e116721ab5.png
   :alt: Final beam distribution when running with full resolution (see inline comments in the input file/script).

   Final beam distribution when running with full resolution (see inline comments in the input file/script).

