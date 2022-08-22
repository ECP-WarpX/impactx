.. _install-dependencies:

Dependencies
============

ImpactX depends on the following popular third party software.
Please see installation instructions below.

- a mature `C++17 <https://en.wikipedia.org/wiki/C%2B%2B17>`__ compiler, e.g., GCC 7, Clang 7, NVCC 11.0, MSVC 19.15 or newer
- `CMake 3.15.0+ <https://cmake.org>`__
- `Git 2.18+ <https://git-scm.com>`__
- `AMReX <https://amrex-codes.github.io>`__: we automatically download and compile a copy
- `WarpX <https://github.com/ECP-WarpX/warpx>`__: we automatically download and compile a copy

Optional dependencies include:

- `MPI 3.0+ <https://www.mpi-forum.org/docs/>`__: for multi-node and/or multi-GPU execution
- `CUDA Toolkit 11.0+ <https://developer.nvidia.com/cuda-downloads>`__: for Nvidia GPU support (see `matching host-compilers <https://gist.github.com/ax3l/9489132>`_)
- `OpenMP 3.1+ <https://www.openmp.org>`__: for threaded CPU execution
- `FFTW3 <http://www.fftw.org>`_: for spectral solver support
- `openPMD-api 0.14.2+ <https://github.com/openPMD/openPMD-api>`__: we automatically download and compile a copy of openPMD-api for openPMD I/O support

  - see `optional I/O backends <https://github.com/openPMD/openPMD-api#dependencies>`__
- `CCache <https://ccache.dev>`__: to speed up rebuilds (needs 3.7.9+ for CUDA)
- `Ninja <https://ninja-build.org>`__: for faster parallel compiles


Install
-------

Pick *one* of the installation methods below to install all dependencies for ImpactX development in a consistent manner.

Spack (macOS/Linux)
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   spack env create impactx-dev
   spack env activate impactx-dev
   spack add adios2        # for openPMD
   spack add ccache
   spack add cmake
   spack add fftw
   spack add hdf5          # for openPMD
   spack add mpi
   spack add pkgconfig     # for fftw
   spack add python
   spack add py-pip
   spack add py-setuptools
   spack add py-wheel

   # OpenMP support on macOS
   [[ $OSTYPE == 'darwin'* ]] && spack add llvm-openmp

   # optional: Linux only
   #spack add cuda

   spack install
   python3 -m pip install matplotlib numpy openpmd-api pandas pytest scipy

In new terminals, re-activate the environment with ``spack env activate impactx-dev`` again.


Brew (macOS/Linux)
^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   brew update
   brew install adios2      # for openPMD
   brew install ccache
   brew install cmake
   brew install fftw
   brew install git
   brew install hdf5-mpi    # for openPMD
   brew install libomp      # for OpenMP
   brew install pkg-config  # for fftw
   brew install open-mpi
   brew install python
   python3 -m pip install matplotlib yt scipy numpy openpmd-api


Conda (Linux/macOS/Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

With MPI (only Linux/macOS):

.. code-block:: bash

   conda create -n impactx-dev -c conda-forge ccache cmake compilers git "openpmd-api=*=mpi_mpich*" python mpich numpy scipy yt "fftw=*=mpi_mpich*" matplotlib mamba ninja numpy pandas pytest scipy
   conda activate impactx-dev

Without MPI:

.. code-block:: bash

   conda create -n impactx-nompi-dev -c conda-forge ccache cmake compilers git openpmd-api python numpy scipy yt fftw matplotlib mamba ninja numpy pandas scipy
   conda activate impactx-nompi-dev

   # compile ImpactX with -DImpactX_MPI=OFF

.. note::

   A general option to deactivate that conda self-activates its base environment.
   This `avoids interference with the system and other package managers <https://collegeville.github.io/CW20/WorkshopResources/WhitePapers/huebl-working-with-multiple-pkg-mgrs.pdf>`__.

   .. code-block:: bash

      conda config --set auto_activate_base false


Apt (Debian/Ubuntu)
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sudo apt update
   sudo apt install build-essential ccache cmake g++ git libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-matplotlib python3-numpy python3-pandas python3-scipy
