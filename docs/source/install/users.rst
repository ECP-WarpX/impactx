.. _install-users:

Users
=====

.. raw:: html

   <style>
   .rst-content .section>img {
       width: 30px;
       margin-bottom: 0;
       margin-top: 0;
       margin-right: 15px;
       margin-left: 15px;
       float: left;
   }
   </style>

Our community is here to help.
Please `report installation problems <https://github.com/ECP-WarpX/impactx/issues/new>`_ in case you should get stuck.

Choose **one** of the installation methods below to get started:

.. only:: html

   .. image:: hpc.svg

HPC Systems
-----------

If want to use ImpactX on a specific high-performance computing (HPC) systems, jump directly to our :ref:`HPC system-specific documentation <install-hpc>`.


.. _install-conda:

.. only:: html

   .. image:: conda.svg

Using the Conda Package
-----------------------

A package for ImpactX is available via the `Conda <https://conda.io>`_ package manager.

.. code-block:: bash

   conda create -n impactx -c conda-forge impactx
   conda activate impactx

Note: the ``impactx`` `conda package <https://anaconda.org/conda-forge/impactx>`__ does not yet provide GPU support.


.. _install-spack:

.. only:: html

   .. image:: spack.svg

Using the Spack Package
-----------------------

.. note::

   Coming soon.


.. _install-pypi:

.. only:: html

   .. image:: pypi.svg

Using the PyPI Package
----------------------

.. note::

   Coming soon.


.. _install-brew:

.. only:: html

   .. image:: brew.svg

Using the Brew Package
----------------------

.. note::

   Coming soon.


.. _install-cmake:

.. only:: html

   .. image:: cmake.svg

From Source with CMake
----------------------

After installing the :ref:`ImpactX dependencies <install-dependencies>`, you can also install ImpactX from source with `CMake <https://cmake.org/>`_:

.. code-block:: bash

   # get the source code
   git clone https://github.com/ECP-WarpX/impactx.git $HOME/src/impactx
   cd $HOME/src/impactx

   # configure
   cmake -S . -B build

   # optional: change configuration
   ccmake build

   # compile
   #   on Windows:          --config Release
   cmake --build build -j 4

   # executables for ImpactX are now in build/bin/

We document the details in the :ref:`developer installation <install-developers>`.

Tips for macOS Users
--------------------

.. tip::

   Before getting started with package managers, please check what you manually installed in ``/usr/local``.
   If you find entries in ``bin/``, ``lib/`` et al. that look like you manually installed MPI, HDF5 or other software in the past, then remove those files first.

   If you find software such as MPI in the same directories that are shown as symbolic links then it is likely you `brew installed <https://brew.sh/>`__ software before.
   If you are trying annother package manager than ``brew``, run `brew unlink ... <https://docs.brew.sh/Tips-N%27-Tricks#quickly-remove-something-from-usrlocal>`__ on such packages first to avoid software incompatibilities.

See also: A. Huebl, `Working With Multiple Package Managers <https://collegeville.github.io/CW20/WorkshopResources/WhitePapers/huebl-working-with-multiple-pkg-mgrs.pdf>`__, `Collegeville Workshop (CW20) <https://collegeville.github.io/CW20/>`_, 2020
