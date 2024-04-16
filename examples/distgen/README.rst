
Generation of beam distributions
=================================

The following examples are tests of beam initialization for distributions of various types.

In each example, we use a 2 GeV electron beam with initial unnormalized rms emittance of 2 nm.

The matched Twiss parameters are the same as those used in the FODO example:

* :math:`\beta_\mathrm{x} = 2.82161941` m
* :math:`\alpha_\mathrm{x} = -1.59050035`
* :math:`\beta_\mathrm{y} = 2.82161941` m
* :math:`\alpha_\mathrm{y} = 1.59050035`

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


.. _examples-gaussian:

A 6d Gaussian distribution
============================

A Gaussian distribution in all 6 phase space variables.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_gaussian.py`` or
* ImpactX **executable** using an input file: ``impactx input_gaussian.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_gaussian.py
          :language: python3
          :caption: You can copy this file from ``examples/distgen/run_gaussian.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_gaussian.in
          :language: ini
          :caption: You can copy this file from ``examples/distgen/input_gaussian.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_gaussian.py``

   .. literalinclude:: analysis_gaussian.py
      :language: python3
      :caption: You can copy this file from ``examples/distgen/analysis_gaussian.py``.



.. _examples-kvdist:

A Kapchinskij-Vladimirskij (K-V) distribution
===============================================

A 4D K-V distribution in the transverse phase space variables ( + a longitudinally uniform distribution in :math:`t` + a Gaussian distribution in :math:`p_t` ).

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_kvdist.py`` or
* ImpactX **executable** using an input file: ``impactx input_kvdist.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_kvdist.py
          :language: python3
          :caption: You can copy this file from ``examples/distgen/run_kvdist.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_kvdist.in
          :language: ini
          :caption: You can copy this file from ``examples/distgen/input_kvdist.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_kvdist.py``

   .. literalinclude:: analysis_kvdist.py
      :language: python3
      :caption: You can copy this file from ``examples/distgen/analysis_kvdist.py``.




.. _examples-kvdist_from_twiss:

A K-V distribution initialized from Twiss functions
======================================================

Identical to the previous example (examples-kvdist), but initialized using Courant-Snyder Twiss functions.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_kvdist_from_twiss.py`` or
* ImpactX **executable** using an input file: ``impactx input_kvdist_from_twiss.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_kvdist_from_twiss.py
          :language: python3
          :caption: You can copy this file from ``examples/distgen/run_kvdist_from_twiss.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_kvdist_from_twiss.in
          :language: ini
          :caption: You can copy this file from ``examples/distgen/input_kvdist_from_twiss.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_kvdist_from_twiss.py``

   .. literalinclude:: analysis_kvdist_from_twiss.py
      :language: python3
      :caption: You can copy this file from ``examples/distgen/analysis_kvdist_from_twiss.py``.




.. _examples-kurth4d:

A 4D Kurth Distribution
============================

A 4D Kurth distribution in the transverse phase space variables ( + a longitudinally uniform distribution in :math:`t` + a Gaussian distribution in :math:`p_t` ).


In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_kurth4d.py`` or
* ImpactX **executable** using an input file: ``impactx input_kurth4d.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_kurth4d.py
          :language: python3
          :caption: You can copy this file from ``examples/distgen/run_kurth4d.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_kurth4d.in
          :language: ini
          :caption: You can copy this file from ``examples/distgen/input_kurth4d.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_kurth4d.py``

   .. literalinclude:: analysis_kurth4d.py
      :language: python3
      :caption: You can copy this file from ``examples/distgen/analysis_kurth4d.py``.




.. _examples-semigaussian:

A Semigaussian distribution
============================

A 6D semigaussian distribution (uniform in position, Gaussian in momentum).

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_semigaussian.py`` or
* ImpactX **executable** using an input file: ``impactx input_semigaussian.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_semigaussian.py
          :language: python3
          :caption: You can copy this file from ``examples/distgen/run_semigaussian.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_semigaussian.in
          :language: ini
          :caption: You can copy this file from ``examples/distgen/input_semigaussian.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_semigaussian.py``

   .. literalinclude:: analysis_semigaussian.py
      :language: python3
      :caption: You can copy this file from ``examples/distgen/analysis_semigaussian.py``.
