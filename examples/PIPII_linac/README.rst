.. _examples-pipii-linac:

The Fermilab PIP-II Linac
===========================================================

A version of the lattice of Fermilab's PIP-II linac, describing acceleration of a 5 mA H- / proton beam from 2.1 MeV to 800 MeV.

This lattice describes propagation from the exit of the RFQ to the entry of the Booster. 

This file was converted from TraceWin input using a Jupyter notebook, whose corresponding script is given as "parse_lattice.py".

The TraceWin lattice parser is based on https://github.com/ChristopherMayes/lume-tracewin .

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as a Python script (``python3 run_pipii_linac.py``) or with an app with an input file (``impactx input_pipii_linac.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_pipii_linac.py
          :language: python3
          :caption: You can copy this file from ``examples/PIPII_linac/run_pipii_linac.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_pipii_linac.in
          :language: ini
          :caption: You can copy this file from ``examples/PIPII_linac/input_pipii_linac.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_pipii_linac.py``

   .. literalinclude:: analysis_pipii_linac.py
      :language: python3
      :caption: You can copy this file from ``examples/PIPII_linac/analysis_pipii_linac.py``.
