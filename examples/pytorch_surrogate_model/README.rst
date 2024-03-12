.. _examples-ml-surrogate:

15 Stage Laser-Plasma Accelerator Surrogate
===========================================

This example models an electron beam accelerated through fifteen stages of laser-plasma accelerators with ideal plasma lenses providing the focusing between stages.
For more details, see:


- Sandberg R T, Lehe R, Mitchell C E, Garten M, Qiang J, Vay J-L and Huebl A.
  **Synthesizing Particle-in-Cell Simulations Through Learning and GPU Computing for Hybrid Particle Accelerator Beamlines**.
  Proc. of Platform for Advanced Scientific Computing (PASC'24), *submitted*, 2024.
- Sandberg R T, Lehe R, Mitchell C E, Garten M, Qiang J, Vay J-L and Huebl A.
  **Hybrid Beamline Element ML-Training for Surrogates in the ImpactX Beam-Dynamics Code**.
  14th International Particle Accelerator Conference (IPAC'23), WEPA101, 2023.
  `DOI:10.18429/JACoW-IPAC2023-WEPA101 <https://doi.org/10.18429/JACoW-IPAC2023-WEPA101>`__

A schematic with more information can be seen in the figure below:

.. figure:: https://private-user-images.githubusercontent.com/10621396/312266448-a51613b6-7ea0-4f5f-b555-66b3db39fefa.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MTAyOTk4ODAsIm5iZiI6MTcxMDI5OTU4MCwicGF0aCI6Ii8xMDYyMTM5Ni8zMTIyNjY0NDgtYTUxNjEzYjYtN2VhMC00ZjVmLWI1NTUtNjZiM2RiMzlmZWZhLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNDAzMTMlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjQwMzEzVDAzMTMwMFomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTZjYmMzNjBlNmUzZDllZjNmODE5YTE0YjU4NzIzNmU5NGQ2YWUxNzEwNWVjOWNlOTY3YzdiNDNkMmI4MmRhMmQmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.lBZTlJNHpCaL8unPe4Hb_DjG25RMlnCUWztv68EqdQs
   :alt: [fig:lpa_schematic] Schematic of the 15 stages of laser-plasma accelerators.

   [fig:lpa_schematic] Schematic of the 15 stages of laser-plasma accelerators.

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

.. figure:: https://user-images.githubusercontent.com/10621396/289976300-6f861d19-5a5c-42eb-9435-9f57bd2010bf.png
   :alt: Evolution of beam moments through 15 stage LPA via neural network surrogates.

   Evolution of electron beam moments through 15 stages of LPAs (via neural network surrogates).

.. figure:: https://user-images.githubusercontent.com/10621396/289956805-49e0a94a-454f-4b48-b448-7ac772edf95a.png
   :alt: Initial phase space projections

   Initial phase space projections going into 15 stage LPA (via neural network surrogates) simulation. Top row: spatial projections, middle row: momentum projections, bottom row: phase spaces.

.. figure:: https://user-images.githubusercontent.com/10621396/289975961-7d389864-9106-4446-8556-b0ea4bb28145.png
   :alt: Final phase space projections after 15 stage LPA (via neural network surrogates) simulation

   Final phase space projections after 15 stage LPA (via neural network surrogates) simulation. Top row: spatial projections, middle row: momentum projections, bottom row: phase spaces.
