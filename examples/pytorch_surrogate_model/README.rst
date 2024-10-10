.. _examples-ml-surrogate:

15 Stage Laser-Plasma Accelerator Surrogate
===========================================

This example models an electron beam accelerated through fifteen stages of laser-plasma accelerators with ideal plasma lenses providing the focusing between stages.
For more details, see:

- Sandberg R T, Lehe R, Mitchell C E, Garten M, Myers A, Qiang J, Vay J-L and Huebl A.
  **Synthesizing Particle-in-Cell Simulations Through Learning and GPU Computing for Hybrid Particle Accelerator Beamlines**.
  Proc. of Platform for Advanced Scientific Computing (PASC'24), *PASC24 Best Paper Award*, 2024.
  `DOI:10.1145/3659914.3659937 <https://doi.org/10.1145/3659914.3659937>`__

- Sandberg R T, Lehe R, Mitchell C E, Garten M, Qiang J, Vay J-L and Huebl A.
  **Hybrid Beamline Element ML-Training for Surrogates in the ImpactX Beam-Dynamics Code**.
  14th International Particle Accelerator Conference (IPAC'23), WEPA101, 2023.
  `DOI:10.18429/JACoW-IPAC2023-WEPA101 <https://doi.org/10.18429/JACoW-IPAC2023-WEPA101>`__

A schematic with more information can be seen in the figure below:

.. figure:: https://gist.githubusercontent.com/RTSandberg/cf3f6193b3da12e7fd815f69789fd6f2/raw/2308e412d7482d55811afa2e9a0c6fa97627fc2f/schema_15_stages.png
   :alt: Schematic of the 15 stages of laser-plasma accelerators.

   Schematic of the 15 stages of laser-plasma accelerators.

The laser-plasma accelerator elements are modeled with neural networks as surrogates.
These networks are trained beforehand.
In this example, pre-trained neural networks are downloaded from `a Zenodo archive <https://zenodo.org/records/10368972>`__ and saved in the ``models`` directory.
For more about how these neural network surrogate models were created,
see this `description of a workflow for training neural networks from WarpX simulation data <https://warpx.readthedocs.io/en/latest/usage/workflows/ml_dataset_training.html>`__.

We use a 1 GeV electron beam with initial normalized rms emittance of 1 mm-mrad.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

Run
---

This example can **only** be run with **Python**:

* **Python** script: ``python3 run_ml_surrogate_15_stage.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: run_ml_surrogate_15_stage.py
   :language: python
   :caption: You can copy this file from ``examples/pytorch_surrogate_model/run_ml_surrogate_15_stage.py``.


This script requires some utility code for using the neural networks that is provided here:

.. dropdown:: Script ``surrogate_model_definitions.py``

   .. literalinclude:: surrogate_model_definitions.py
      :language: python
      :caption: You can copy this file from ``examples/pytorch_surrogate_model/surrogate_model_definitions.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analyze_ml_surrogate_15_stage.py``

   .. literalinclude:: analyze_ml_surrogate_15_stage.py
      :language: python
      :caption: You can copy this file from ``examples/pytorch_surrogate_model/analyze_ml_surrogate_15_stage.py``.

Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``visualize_ml_surrogate_15_stage.py``

   .. literalinclude:: visualize_ml_surrogate_15_stage.py
      :language: python
      :caption: You can copy this file from ``examples/pytorch_surrogate_model/visualize_ml_surrogate_15_stage.py``.

.. figure:: https://gist.githubusercontent.com/RTSandberg/cf3f6193b3da12e7fd815f69789fd6f2/raw/2308e412d7482d55811afa2e9a0c6fa97627fc2f/lpa_ml_surrogate_moments.png
   :alt: Evolution of beam moments through 15 stage LPA via neural network surrogates.

   Evolution of electron beam moments through 15 stages of LPAs (via neural network surrogates).

.. figure:: https://gist.githubusercontent.com/RTSandberg/cf3f6193b3da12e7fd815f69789fd6f2/raw/2308e412d7482d55811afa2e9a0c6fa97627fc2f/initial_phase_spaces.png
   :alt: Initial phase space projections

   Initial phase space projections going into 15 stage LPA (via neural network surrogates) simulation. Top row: spatial projections, middle row: momentum projections, bottom row: phase spaces.

.. figure:: https://gist.githubusercontent.com/RTSandberg/cf3f6193b3da12e7fd815f69789fd6f2/raw/2308e412d7482d55811afa2e9a0c6fa97627fc2f/stage_14_phase_spaces.png
   :alt: Final phase space projections after 15 stage LPA (via neural network surrogates) simulation

   Final phase space projections after 15 stage LPA (via neural network surrogates) simulation. Top row: spatial projections, middle row: momentum projections, bottom row: phase spaces.
