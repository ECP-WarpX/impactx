.. _building-lumi:

LUMI (CSC)
==========

The `LUMI cluster <https://www.lumi-supercomputer.eu>`__ is located at CSC (Finland).
LUMI has multiple partitions.

On **LUMI-G**, each node contains `four (4) AMD MI250X GPUs <https://docs.lumi-supercomputer.eu/hardware/lumig/>`__, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

On **LUMI-C**, each node contains `two (2) AMD EPYC 7763 CPUs <https://docs.lumi-supercomputer.eu/hardware/lumic/>`__.


Introduction
------------

If you are new to this system, **please see the following resources**:

* `Lumi user guide <https://docs.lumi-supercomputer.eu>`__

  * `Project Maintainance <https://my.lumi-supercomputer.eu>`__ and `SSH Key management <https://mms.myaccessid.org>`__
  * `Quotas and projects <https://docs.lumi-supercomputer.eu/runjobs/lumi_env/dailymanagement/>`__
* Batch system: `Slurm <https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/slurm-quickstart/>`__
* `Data analytics and visualization <https://docs.lumi-supercomputer.eu/hardware/lumid/>`__
* `Production directories <https://docs.lumi-supercomputer.eu/storage/>`__:

  * ``$HOME``: single user, intended to store user configuration files and personal data (20GB default quota)
  * ``/project/$proj``: shared with all members of a project, purged at the end of a project (50 GB default quota)
  * ``/scratch/$proj``: temporary storage, main storage to be used for disk I/O needs when running simulations on LUMI, purged every 90 days (50TB default quota)


.. _building-lumi-preparation:

Preparation
-----------

Use the following commands to download the ImpactX source code:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/impactx.git $HOME/src/impactx

On LUMI, you can run either with fast MI250X GPUs (LUMI-G, recommended) or CPU nodes (LUMI-C).

.. tab-set::

   .. tab-item:: MI250X GPUs

      We use system software modules, add environment hints and further dependencies via the file ``$HOME/lumi_gpu_impactx.profile``.
      Create it now:

      .. code-block:: bash

         cp $HOME/src/impactx/docs/source/install/hpc/lumi-csc/lumi_gpu_impactx.profile.example $HOME/lumi_gpu_impactx.profile

      .. dropdown:: Script Details
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: lumi-csc/lumi_gpu_impactx.profile.example
            :language: bash

      Edit the 2nd line of this script, which sets the ``export proj="project_..."`` variable using a text editor
      such as ``nano``, ``emacs``, or ``vim`` (all available by default on LUMI login nodes).
      You can find out your project name by running ``lumi-ldap-userinfo`` on LUMI.
      For example, if you are member of the project ``project_465000962``, then run ``nano $HOME/lumi_gpu_impactx.profile`` and edit line 2 to read:

      .. code-block:: bash

         export proj="project_465000962"

      Exit the ``nano`` editor with ``Ctrl`` + ``O`` (save) and then ``Ctrl`` + ``X`` (exit).

      .. important::

         Now, and as the first step on future logins to LUMI, activate these environment settings:

         .. code-block:: bash

            source $HOME/lumi_gpu_impactx.profile

      Finally, since LUMI does not yet provide software modules for some of our dependencies, install them once:

      .. code-block:: bash

         bash $HOME/src/impactx/docs/source/install/hpc/lumi-csc/install_gpu_dependencies.sh
         source $HOME/sw/lumi/gpu/venvs/impactx-gpu-lumi/bin/activate

      .. dropdown:: Script Details
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: lumi-csc/install_gpu_dependencies.sh
            :language: bash

   .. tab-item:: CPU Nodes

      We use system software modules, add environment hints and further dependencies via the file ``$HOME/lumi_cpu_impactx.profile``.
      Create it now:

      .. code-block:: bash

         cp $HOME/src/impactx/docs/source/install/hpc/lumi-csc/lumi_cpu_impactx.profile.example $HOME/lumi_cpu_impactx.profile

      .. dropdown:: Script Details
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: lumi-csc/lumi_cpu_impactx.profile.example
            :language: bash

      Edit the 2nd line of this script, which sets the ``export proj="project_..."`` variable using a text editor
      such as ``nano``, ``emacs``, or ``vim`` (all available by default on LUMI login nodes).
      You can find out your project name by running ``lumi-ldap-userinfo`` on LUMI.
      For example, if you are member of the project ``project_465000962``, then run ``nano $HOME/lumi_cpu_impactx.profile`` and edit line 2 to read:

      .. code-block:: bash

         export proj="project_465000962"

      Exit the ``nano`` editor with ``Ctrl`` + ``O`` (save) and then ``Ctrl`` + ``X`` (exit).

      .. important::

         Now, and as the first step on future logins to LUMI, activate these environment settings:

         .. code-block:: bash

            source $HOME/lumi_cpu_impactx.profile

      Finally, since LUMI does not yet provide software modules for some of our dependencies, install them once:

      .. code-block:: bash

         bash $HOME/src/impactx/docs/source/install/hpc/lumi-csc/install_cpu_dependencies.sh
         source $HOME/sw/lumi/cpu/venvs/impactx-cpu-lumi/bin/activate

      .. dropdown:: Script Details
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: lumi-csc/install_cpu_dependencies.sh
            :language: bash


.. _building-lumi-compilation:

Compilation
-----------

Use the following :ref:`cmake commands <building-cmake>` to compile the application executable:

.. tab-set::

   .. tab-item:: MI250X GPUs

      .. code-block:: bash

         cd $HOME/src/impactx

         cmake --fresh -S . -B build_lumi_gpu -DImpactX_COMPUTE=HIP -DImpactX_FFT=ON
         cmake --build build_lumi_gpu -j 16

      The ImpactX application executables are now in ``$HOME/src/impactx/build_lumi_gpu/bin/``.
      Additionally, the following commands will install ImpactX as a Python module:

      .. code-block:: bash

         cmake --fresh -S . -B build_lumi_gpu_py -DImpactX_COMPUTE=HIP -DImpactX_FFT=ON -DImpactX_PYTHON=ON
         cmake --build build_lumi_gpu_py -j 16 --target pip_install

   .. tab-item:: CPU Nodes

      .. code-block:: bash

         cd $HOME/src/impactx

         cmake --fresh -S . -B build_lumi_cpu -DImpactX_FFT=ON
         cmake --build build_lumi_cpu -j 16

      The ImpactX application executables are now in ``$HOME/src/impactx/build_lumi_cpu/bin/``.
      Additionally, the following commands will install ImpactX as a Python module:

      .. code-block:: bash

         cmake --fresh -S . -B build_lumi_cpu_py -DImpactX_FFT=ON -DImpactX_PYTHON=ON
         cmake --build build_lumi_cpu_py -j 16 --target pip_install

Now, you can :ref:`submit LUMI compute jobs <running-cpp-lumi>` for ImpactX :ref:`Python scripts <usage-picmi>` (:ref:`example scripts <usage-examples>`).
Or, you can use the ImpactX executables to submit LUMI jobs (:ref:`example inputs <usage-examples>`).
For executables, you can reference their location in your :ref:`job script <running-cpp-lumi>` or copy them to a location in ``/project/$proj`` or ``/scratch/$proj``.


.. _building-lumi-update:

Update ImpactX & Dependencies
-----------------------------

If you already installed ImpactX in the past and want to update it, start by getting the latest source code:

.. code-block:: bash

   cd $HOME/src/impactx

   # read the output of this command - does it look ok?
   git status

   # get the latest ImpactX source code
   git fetch
   git pull

   # read the output of these commands - do they look ok?
   git status
   git log     # press q to exit

And, if needed,

- :ref:`update the lumi_gpu_impactx.profile and lumi_cpu_impactx.profile files <building-lumi-preparation>`,
- log out and into the system, activate the now updated environment profile as usual,
- :ref:`execute the dependency install scripts <building-lumi-preparation>`.

As the last step, :ref:`recompile ImpactX <building-lumi-compilation>`.


.. _running-cpp-lumi:

Running
-------

On LUMI, compute jobs are run via the `Slurm <https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/slurm-quickstart/>`__ resource manager.
There are various `Slurm priority queues <https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions/>`__ besides the defaults documented below.

For interactive runs on a single node, simply use the aliases ``getNode`` or ``runNode ...`` that are defined by the ``lumi_*_impactx.profile``.
For runs on multiple nodes, use and adjust the templates below.

.. tab-set::

   .. tab-item:: MI250X GPUs (2x64 GB)

      The **LUMI-G** (GPU) partition on the supercomputer LUMI at CSC has up to `2978 nodes <https://docs.lumi-supercomputer.eu/hardware/lumig/>`__, each with 8 Graphics Compute Dies (GCDs).
      ImpactX runs one MPI rank per Graphics Compute Die.

      The batch script below can be used to run an ImpactX simulation on multiple nodes (change ``-N`` accordingly).
      Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<project id>`` or the concete inputs file.
      Copy the executable or point to it via ``EXE`` and adjust the path for the ``INPUTS`` variable accordingly.

      .. literalinclude:: lumi-csc/lumi_gpu.sbatch
         :language: bash
         :caption: You can copy this file from ``$HOME/src/impactx/docs/source/install/hpc/lumi-csc/lumi_gpu.sbatch``.

      To run a simulation, copy the lines above to a file ``lumi_gpu.sbatch`` and run

      .. code-block:: bash

         sbatch lumi_gpu.sbatch

      to submit the job.

   .. tab-item:: CPU Nodes

      The **LUMI-C** (CPU) partition on the supercomputer LUMI at CSC has up to `2048 nodes <https://docs.lumi-supercomputer.eu/hardware/lumic/>`__.
      ImpactX runs 16 MPI ranks per node, each with 8 OpenMP threads.

      The batch script below can be used to run an ImpactX simulation on multiple nodes (change ``-N`` accordingly).
      Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<project id>`` or the concete inputs file.
      Copy the executable or point to it via ``EXE`` and adjust the path for the ``INPUTS`` variable accordingly.

      .. literalinclude:: lumi-csc/lumi_cpu.sbatch
         :language: bash
         :caption: You can copy this file from ``$HOME/src/impactx/docs/source/install/hpc/lumi-csc/lumi_cpu.sbatch``.

      To run a simulation, copy the lines above to a file ``lumi_cpu.sbatch`` and run

      .. code-block:: bash

         sbatch lumi_cpu.sbatch

      to submit the job.


.. _post-processing-lumi:

Post-Processing
---------------

LUMI provides a `Jupyter Lab service <https://docs.lumi-supercomputer.eu/runjobs/webui/jupyter/>`__ that can be used for interactive post-processing.

You can reuse your CPU ImpactX virtual environment for post-processing in Jupyter Lab:
* Python: ``Cray-python (3.10.10)``
* Virtual environment path: ``${HOME}/sw/lumi/cpu/venvs/impactx-cpu-lumi/`` (replace ``${HOME}`` with the output of ``echo ${HOME}``)


.. _known-lumi-issues:

Known System Issues
-------------------

.. warning::

   December 12th, 2022:
   There is a caching bug in libFabric that causes ImpactX simulations to occasionally hang on LUMI on more than 1 node.

   As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

   .. code-block:: bash

      #export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
      # or, less invasive:
      export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

.. warning::

   January, 2023:
   We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.

   June, 2023:
   Although a fix was planned for ROCm 5.5, we still see the same issue in this release and continue to exchange with AMD and HPE on the issue.

   Stay with the ROCm 5.2 module to avoid a 2x slowdown.

.. warning::

   May 2023:
   rocFFT in ROCm 5.1-5.3 tries to `write to a cache <https://rocfft.readthedocs.io/en/latest/#runtime-compilation>`__ in the home area by default.
   This does not scale, disable it via:

   .. code-block:: bash

      export ROCFFT_RTC_CACHE_PATH=/dev/null
