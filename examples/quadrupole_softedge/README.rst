.. _examples-quadrupole-softedge:

Soft-Edge Quadrupole
===================

Proton beam propagating through a 6 m region containing a soft-edge
quadrupole.

The quadrupole model used is the default thin-shell model described in (to modify):
P. Granum et al, "Efficient calculations of magnetic fields of solenoids for simulations,"
NIMA 1034, 166706 (2022)
`DOI:10.1016/j.nima.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1 micron
in the horizontal plane, and 2 micron in the vertical plane.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as a Python script (``python3 run_quadrupole_softedge.py``) or with an app with an input file (``impactx input_quadrupole_softedge.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_quadrupole_softedge.py
          :language: python3
          :caption: You can copy this file from ``examples/quadrupole_softedge/run_quadrupole_softedge.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_quadrupole_softedge.in
          :language: ini
          :caption: You can copy this file from ``examples/quadrupole_softedge/input_quadrupole_softedge.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_quadrupole_softedge.py``

   .. literalinclude:: analysis_quadrupole_softedge.py
      :language: python3
      :caption: You can copy this file from ``examples/quadrupole_softedge/analysis_quadrupole_softedge.py``.
