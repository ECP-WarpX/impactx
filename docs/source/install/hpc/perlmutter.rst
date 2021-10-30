.. _building-perlmutter:

Perlmutter (NERSC)
==================

.. warning::

   Perlmutter is still in acceptance testing.
   This page documents our internal testing workflow only.

The `Perlmutter cluster <https://docs.nersc.gov/systems/perlmutter/>`_ is located at NERSC.

If you are new to this system, please see the following resources:

* `NERSC user guide <https://docs.nersc.gov/>`__
* Batch system: `Slurm <https://docs.nersc.gov/systems/perlmutter/#running-jobs>`__
* `Jupyter service <https://docs.nersc.gov/services/jupyter/>`__
* `Production directories <https://docs.nersc.gov/filesystems/perlmutter-scratch/>`__:

  * ``$PSCRATCH``: per-user production directory (<TBD>TB)
  * ``/global/cscratch1/sd/m3239``: shared production directory for users in the project ``m3239`` (50TB)
  * ``/global/cfs/cdirs/m3239/``: community file system for users in the project ``m3239`` (100TB)


Installation
------------

Use the following commands to download the ImpactX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/impactx.git $HOME/src/impactx

We use the following modules and environments on the system (``$HOME/perlmutter_impactx.profile``).

.. code-block:: bash

   # please set your project account
   export proj=<yourProject>  # LBNL/AMP: m3906_g

   # required dependencies
   module load cmake/git-20210830  # 3.22-dev
   module swap PrgEnv-nvidia PrgEnv-gnu
   module swap gcc gcc/9.3.0
   module load cuda

   # optional: just an additional text editor
   # module load nano  # TODO: request from support

   # optional: Python, ...
   # TODO

   # optional: an alias to request an interactive node for two hours
   function getNode() {
       salloc -N 1 --ntasks-per-node=4 -t 2:00:00 -C gpu -c 32 -G 4 -A $proj
   }

   # GPU-aware MPI
   export MPICH_GPU_SUPPORT_ENABLED=1

   # optimize CUDA compilation for A100
   export AMREX_CUDA_ARCH=8.0

   # compiler environment hints
   export CC=$(which gcc)
   export CXX=$(which g++)
   export FC=$(which gfortran)
   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which g++)


We recommend to store the above lines in a file, such as ``$HOME/perlmutter_impactx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/perlmutter_impactx.profile

Then, ``cd`` into the directory ``$HOME/src/impactx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/impactx
   rm -rf build

   cmake -S . -B build -DImpactX_OPENPMD=ON -DImpactX_COMPUTE=CUDA
   cmake --build build -j 32

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.


.. _running-cpp-perlmutter:

Running
-------

.. _running-cpp-perlmutter-A100-GPUs:

A100 GPUs (40 GB)
^^^^^^^^^^^^^^^^^

The batch script below can be used to run a ImpactX simulation on multiple nodes (change ``-N`` accordingly) on the supercomputer Perlmutter at NERSC.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
Note that we run one MPI rank per GPU.


.. literalinclude:: ../../../../etc/impactx/perlmutter-nersc/batch_perlmutter.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_perlmutter.sh`` and run

.. code-block:: bash

   sbatch batch_perlmutter.sh

to submit the job.


.. _post-processing-perlmutter:

Post-Processing
---------------

For post-processing, most users use Python via NERSC's `Jupyter service <https://jupyter.nersc.gov>`__ (`Docs <https://docs.nersc.gov/services/jupyter/>`__).

Please follow the same guidance as for :ref:`NERSC Cori post-processing <post-processing-cori>`.

The Perlmutter ``$PSCRATCH`` filesystem is currently not yet available on Jupyter.
Thus, store or copy your data to Cori's ``$SCRATCH`` or use the Community FileSystem (CFS) for now.
