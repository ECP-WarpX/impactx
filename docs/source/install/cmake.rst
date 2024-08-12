.. _install-developers:
.. _building-cmake:
.. _building-cmake-intro:

Developers
==========

`CMake <https://cmake.org>`_ is our primary build system.
If you are new to CMake, `this short tutorial <https://hsf-training.github.io/hsf-training-cmake-webpage/>`_ from the HEP Software foundation is the perfect place to get started.
If you just want to use CMake to build the project, jump into sections `1. Introduction <https://hsf-training.github.io/hsf-training-cmake-webpage/01-intro/index.html>`__, `2. Building with CMake <https://hsf-training.github.io/hsf-training-cmake-webpage/02-building/index.html>`__ and `9. Finding Packages <https://hsf-training.github.io/hsf-training-cmake-webpage/09-findingpackages/index.html>`__.

Dependencies
------------

Before you start, you will need a copy of the ImpactX source code:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/impactx.git $HOME/src/impactx
   cd $HOME/src/impact

ImpactX depends on popular third party software.

* On your development machine, :ref:`follow the instructions here <install-dependencies>`.
* If you are on an HPC machine, :ref:`follow the instructions here <install-hpc>`.

.. toctree::
   :hidden:

   dependencies

.. note::

   Preparation: make sure you work with up-to-date Python tooling.

   .. code-block:: bash

      python3 -m pip install -U pip
      python3 -m pip install -U build packaging setuptools wheel pytest
      python3 -m pip install -U -r examples/requirements.txt


Compile
-------

From the base of the ImpactX source directory, execute:

.. code-block:: bash

   # find dependencies & configure
   #   see additional options below, e.g.
   #                   -DCMAKE_INSTALL_PREFIX=$HOME/sw/impactX
   cmake -S . -B build -DImpactX_PYTHON=ON

   # compile, here we use four threads
   cmake --build build -j 4

That's all! ImpactX binaries are now in ``build/bin/``.
Most people execute these binaries directly or copy them out.

If you want to install the executables in a programmatic way, run this:

.. code-block:: bash

   # for default install paths, you will need administrator rights, e.g. with sudo:
   #   this installs the application
   cmake --build build --target install

   #   this installs the Python bindings via "python3 -m pip install ..."
   cmake --build build --target pip_install -j 4

You can inspect and modify build options after running ``cmake -S . -B build`` with either

.. code-block:: bash

   ccmake build

or by adding arguments with ``-D<OPTION>=<VALUE>`` to the first CMake call, e.g.:

.. code-block:: bash

   cmake -S . -B build -DImpactX_PYTHON=ON -DImpactX_COMPUTE=CUDA

**That's it!**
You can now :ref:`run a first example <usage-examples>`.

Developers could now change the ImpactX source code and then call the install lines again to refresh the installation.

.. tip::

   If you do *not* develop with :ref:`a user-level package manager <install-dependencies>`, e.g., because you rely on a HPC system's environment modules, then consider to set up a virtual environment via `Python venv <https://docs.python.org/3/library/venv.html>`__.
   Otherwise, without a virtual environment, you likely need to add the CMake option ``-DPY_PIP_INSTALL_OPTIONS="--user"``.


Build Options
-------------

=============================== ============================================ ===========================================================
CMake Option                    Default & Values                             Description
=============================== ============================================ ===========================================================
``BUILD_TESTING``               **ON**/OFF                                   Build tests
``CMAKE_BUILD_TYPE``            RelWithDebInfo/**Release**/Debug             Type of build, symbols & optimizations
``CMAKE_INSTALL_PREFIX``        system-dependent path                        Install path prefix
``CMAKE_VERBOSE_MAKEFILE``      ON/**OFF**                                   Print all compiler commands to the terminal during build
``ImpactX_APP``                 **ON**/OFF                                   Build the ImpactX executable application
``ImpactX_COMPUTE``             NOACC/**OMP**/CUDA/SYCL/HIP                  On-node, accelerated computing backend
``ImpactX_FFT``                 ON/**OFF**                                   FFT-based solvers (IGF Space Charge, CSR, ...)
``ImpactX_IPO``                 ON/**OFF**                                   Compile ImpactX with interprocedural optimization (aka LTO)
``ImpactX_MPI``                 **ON**/OFF                                   Multi-node support (message-passing)
``ImpactX_MPI_THREAD_MULTIPLE`` **ON**/OFF                                   MPI thread-multiple support, i.e. for ``async_io``
``ImpactX_OPENPMD``             **ON**/OFF                                   openPMD I/O (HDF5, ADIOS)
``ImpactX_PRECISION``           SINGLE/**DOUBLE**                            Floating point precision (single/double)
``ImpactX_PYTHON``              ON/**OFF**                                   Python bindings
``Python_EXECUTABLE``           (newest found)                               Path to Python executable
``PY_PIP_OPTIONS``              ``-v``                                       Additional options for ``pip``, e.g., ``-vvv``
``PY_PIP_INSTALL_OPTIONS``                                                   Additional options for ``pip install``, e.g., ``--user``
=============================== ============================================ ===========================================================

ImpactX can be configured in further detail with options from AMReX, which are `documented in the AMReX manual <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`_.

**Developers** might be interested in additional options that control dependencies of ImpactX.
By default, the most important dependencies of ImpactX are automatically downloaded for convenience:

============================= ============================================== ===========================================================
CMake Option                  Default & Values                               Description
============================= ============================================== ===========================================================
``BUILD_SHARED_LIBS``         ON/**OFF**                                     Build shared libraries for dependencies
``ImpactX_CCACHE``            **ON**/OFF                                     Search and use CCache to speed up rebuilds.
``ImpactX_ablastr_src``       *None*                                         Path to ABLASTR source directory (preferred if set)
``ImpactX_ablastr_repo``      ``https://github.com/ECP-WarpX/WarpX.git``     Repository URI to pull and build ABLASTR from
``ImpactX_ablastr_branch``    *we set and maintain a compatible commit*      Repository branch for ``ImpactX_ablastr_repo``
``ImpactX_ablastr_internal``  **ON**/OFF                                     Needs a pre-installed ABLASTR library if set to ``OFF``
``ImpactX_amrex_src``         *None*                                         Path to AMReX source directory (preferred if set)
``ImpactX_amrex_repo``        ``https://github.com/AMReX-Codes/amrex.git``   Repository URI to pull and build AMReX from
``ImpactX_amrex_branch``      *we set and maintain a compatible commit*      Repository branch for ``ImpactX_amrex_repo``
``ImpactX_amrex_internal``    **ON**/OFF                                     Needs a pre-installed AMReX library if set to ``OFF``
``ImpactX_openpmd_src``       *None*                                         Path to openPMD-api source directory (preferred if set)
``ImpactX_openpmd_repo``      ``https://github.com/openPMD/openPMD-api.git`` Repository URI to pull and build openPMD-api from
``ImpactX_openpmd_branch``    *we set and maintain a compatible commit*                                     Repository branch for ``ImpactX_openpmd_repo``
``ImpactX_openpmd_internal``  **ON**/OFF                                     Needs a pre-installed openPMD-api library if set to ``OFF``
``ImpactX_pyamrex_src``       *None*                                         Path to AMReX source directory (preferred if set)
``ImpactX_pyamrex_repo``      ``https://github.com/AMReX-Codes/pyamrex.git`` Repository URI to pull and build pyAMReX from
``ImpactX_pyamrex_branch``    *we set and maintain a compatible commit*      Repository branch for ``ImpactX_pyamrex_repo``
``ImpactX_pyamrex_internal``  **ON**/OFF                                     Needs a pre-installed pyAMReX module if set to ``OFF``
``ImpactX_PYTHON_IPO``        **ON**/OFF                                     Build Python w/ interprocedural/link optimization (IPO/LTO)
``ImpactX_pybind11_src``      *None*                                         Path to pybind11 source directory (preferred if set)
``ImpactX_pybind11_repo``     ``https://github.com/pybind/pybind11.git``     Repository URI to pull and build pybind11 from
``ImpactX_pybind11_branch``   *we set and maintain a compatible commit*      Repository branch for ``ImpactX_pybind11_repo``
``ImpactX_pybind11_internal`` **ON**/OFF                                     Needs a pre-installed pybind11 library if set to ``OFF``
============================= ============================================== ===========================================================

For example, one can also build against a local AMReX copy.
Assuming AMReX' source is located in ``$HOME/src/amrex``, add the ``cmake`` argument ``-DImpactX_amrex_src=$HOME/src/amrex``.
Relative paths are also supported, e.g. ``-DImpactX_amrex_src=../amrex``.

Or build against an AMReX feature branch of a colleague.
Assuming your colleague pushed AMReX to ``https://github.com/WeiqunZhang/amrex/`` in a branch ``new-feature`` then pass to ``cmake`` the arguments: ``-DImpactX_amrex_repo=https://github.com/WeiqunZhang/amrex.git -DImpactX_amrex_branch=new-feature``.

If you want to develop against local versions of ABLASTR (from WarpX) and AMReX at the same time, pass for instance ``-DImpactX_ablastr_src=$HOME/src/warpx -DImpactX_amrex_src=$HOME/src/amrex``.

You can speed up the install further if you pre-install these dependencies, e.g. with a package manager.
Set ``-DImpactX_<dependency-name>_internal=OFF`` and add installation prefix of the dependency to the environment variable `CMAKE_PREFIX_PATH <https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html>`__.
Please see the :ref:`introduction to CMake <building-cmake-intro>` if this sounds new to you.

If you re-compile often, consider installing the `Ninja <https://github.com/ninja-build/ninja/wiki/Pre-built-Ninja-packages>`__ build system.
Pass ``-G Ninja`` to the CMake configuration call to speed up parallel compiles.


Configure Your Compiler
-----------------------

If you don't want to use your default compiler, you can set the following environment variables.
For example, using a Clang/LLVM:

.. code-block:: bash

   export CC=$(which clang)
   export CXX=$(which clang++)

If you also want to select a CUDA compiler:

.. code-block:: bash

   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which clang++)

.. note::

   Please clean your build directory with ``rm -rf build/`` after changing the compiler.
   Now call ``cmake -S . -B build`` (+ further options) again to re-initialize the build configuration.


Run
---

The ImpactX Python bindings, which provide the imports ``impactx`` and ``amrex`` (from `pyAMReX <https://github.com/AMReX-Codes/pyamrex>`__), are automatically packaged and installed when calling the ``pip_install`` *CMake target*.

An executable ImpactX application binary with the current compile-time options encoded in its file name will be created in ``build/bin/``.
Additionally, a `symbolic link <https://en.wikipedia.org/wiki/Symbolic_link>`_ named ``impactx`` can be found in that directory, which points to the last built ImpactX executable.
