.. _examples-kicker:

Test of a Transverse Kicker
===========================

This test applies two transverse momentum kicks, first in the horizontal direction (2 mrad) and then in the vertical direction (3 mrad).

We use a 2 GeV electron beam.

The second beam moments should be unchanged, but the first beam moments corresponding to :math:`p_x` and :math:`p_y` should change according to the size of the kick.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_kicker.py`` or
* ImpactX **executable** using an input file: ``impactx input_kicker.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_kicker.py
          :language: python3
          :caption: You can copy this file from ``examples/kicker/run_kicker.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_kicker.in
          :language: ini
          :caption: You can copy this file from ``examples/kicker/input_kicker.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_kicker.py``

   .. literalinclude:: analysis_kicker.py
      :language: python3
      :caption: You can copy this file from ``examples/kicker/analysis_kicker.py``.
