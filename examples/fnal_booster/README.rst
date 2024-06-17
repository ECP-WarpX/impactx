.. _examples-booster:

FNAL Booster Model
===================

This example is based on a thin-kick model of the Fermilab Booster ring originally created in MAD-X, 
and contributed by F. Schmidt (CERN).  The lattice file is expressed in SXF format, and it is parsed and 
executed in ImpactX using a Python script.  A complete set of input and output files can be found on Zenodo:

https://doi.org/10.5281/zenodo.11645618

We use a proton beam with total energy of 1.338272088 GeV and an initial normalized rms emittance of 1.6 microns.

The second moments of the particle distribution after one turn should coincide with the initial second moments of the particle distribution, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 fnal_booster_no-sc.py`` or

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: fnal_booster_no-sc.py
          :language: python3
          :caption: You can copy this file from ``examples/fnal_booster/fnal_booster_no-sc.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_booster.py``

   .. literalinclude:: analysis_booster.py
      :language: python3
      :caption: You can copy this file from ``examples/fnal_booster/analysis_booster.py``.
