.. _examples-ml-surrogate:

9 Stage Laser-Plasma Accelerator Surrogate
==========================================

This example models an electron beam accelerated through nine stages of laser-plasma accelerators with ideal plasma lenses providing the focusing between stages.
For more details, see:


- Sandberg R T, Lehe R, Mitchell C E, Garten M, Qiang J, Vay J-L and Huebl A.
  **Synthesizing Particle-in-Cell Simulations Through Learning and GPU Computing for Hybrid Particle Accelerator Beamlines**.
  Proc. of Platform for Advanced Scientific Computing (PASC'24), *submitted*, 2024.
- Sandberg R T, Lehe R, Mitchell C E, Garten M, Qiang J, Vay J-L and Huebl A.
  **Hybrid Beamline Element ML-Training for Surrogates in the ImpactX Beam-Dynamics Code**.
  14th International Particle Accelerator Conference (IPAC'23), WEPA101, 2023.
  `DOI:10.18429/JACoW-IPAC2023-WEPA101 <https://doi.org/10.18429/JACoW-IPAC2023-WEPA101>`__

A schematic with more information can be seen in the figure below:

.. figure:: https://user-images.githubusercontent.com/10621396/289956389-c7463b99-fb56-490a-8511-22c43f45cdf8.png
   :alt: [fig:lpa_schematic] Schematic of the 9 stages of laser-plasma accelerators.

   [fig:lpa_schematic] Schematic of the 9 stages of laser-plasma accelerators.

The laser-plasma accelerator elements are modeled with neural network surrogates,
previously trained and included in ``models``.
The neural networks require normalized input data; the normalizations can be found in ``datasets``.

We use a 1 GeV electron beam with initial normalized rms emittance of 1 mm-mrad.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

Run
---

This example can be **only** be run with **Python**:

* **Python** script: ``python3 run_ml_surrogate.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: run_ml_surrogate.py
   :language: python
   :caption: You can copy this file from ``examples/pytorch_surrogate_model/run_ml_surrogate.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analyze_ml_surrogate.py``

   .. literalinclude:: analyze_ml_surrogate.py
      :language: python
      :caption: You can copy this file from ``examples/pytorch_surrogate_model/analyze_ml_surrogate.py``.

Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``visualize_ml_surrogate.py``

   .. literalinclude:: visualize_ml_surrogate.py
      :language: python
      :caption: You can copy this file from ``examples/pytorch_surrogate_model/visualize_ml_surrogate.py``.

.. figure:: https://user-images.githubusercontent.com/10621396/289976300-6f861d19-5a5c-42eb-9435-9f57bd2010bf.png
   :alt: Evolution of beam moments through 9 stage LPA via neural network surrogates.

   Evolution of electron beam moments through 9 stages of LPAs (via neural network surrogates).

.. figure:: https://user-images.githubusercontent.com/10621396/289956805-49e0a94a-454f-4b48-b448-7ac772edf95a.png
   :alt: Initial phase space projections

   Initial phase space projections going into 9 stage LPA (via neural network surrogates) simulation. Top row: spatial projections, middle row: momentum projections, bottom row: phase spaces.

.. figure:: https://user-images.githubusercontent.com/10621396/289975961-7d389864-9106-4446-8556-b0ea4bb28145.png
   :alt: Final phase space projections after 9 stage LPA (via neural network surrogates) simulation

   Final phase space projections after 9 stage LPA (via neural network surrogates) simulation. Top row: spatial projections, middle row: momentum projections, bottom row: phase spaces.