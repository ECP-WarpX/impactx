.. _examples-chicane-csr:

Chicane with CSR
================

This is the :ref:`Berlin-Zeuthen magnetic bunch compression chicane <examples-chicane>` example, but this time with coherent synchrotron radiation (CSR) modelled in the bending magnets.

`All parameters can be found online <https://www.desy.de/csr/>`__.
A 5 GeV electron bunch with normalized transverse rms emittance of 1 um undergoes longitudinal compression by a factor of 10 in a standard 4-bend chicane.
The rms pulse length should decrease by the compression factor (10).

An ultrarelativistic model of steady-state CSR in the bending dipoles is included, resulting in a horizontal emittance growth of 19%.  Note that this value is smaller than the horizontal emittance growth of 57% that is obtained when transient (bend-entry and -exit) effects are included.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_chicane_csr.py`` or
* ImpactX **executable** using an input file: ``impactx input_chicane_csr.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_chicane_csr.py
          :language: python3
          :caption: You can copy this file from ``examples/chicane/run_chicane_csr.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_chicane_csr.in
          :language: ini
          :caption: You can copy this file from ``examples/chicane/input_chicane_csr.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_chicane_csr.py``

   .. literalinclude:: analysis_chicane_csr.py
      :language: python3
      :caption: You can copy this file from ``examples/chicane/analysis_chicane_csr.py``.
