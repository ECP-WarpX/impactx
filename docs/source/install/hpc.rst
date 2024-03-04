.. _install-hpc:

HPC
===

On selected high-performance computing (HPC) systems, ImpactX has documented or even pre-build installation routines.
Follow the guide here instead of the generic installation routines for optimal stability and best performance.

impactx.profile
---------------

Use a ``impactx.profile`` file to set up your software environment without colliding with other software.
Ideally, store that file directly in your ``$HOME/`` and source it after connecting to the machine:

.. code-block:: bash

   source $HOME/impactx.profile

We list example ``impactx.profile`` files below, which can be used to set up ImpactX on various HPC systems.

HPC Systems
-----------

.. toctree::
   :maxdepth: 1

   hpc/perlmutter

.. tip::

   Your HPC system is not in the list?
   Our instructions are nearly identical to `installing WarpX, documented here <https://warpx.readthedocs.io/en/latest/install/hpc.html#hpc-machines>`__.

   Also, please do not hesitate to `open an issue <https://github.com/ECP-WarpX/impactx/issues>`__ and together we can document how to run on your preferred system!
