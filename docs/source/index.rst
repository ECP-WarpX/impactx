:orphan:

ImpactX
-------

ImpactX enables high-performance modeling of beam dynamics in particle accelerators with collective effects.

This is the next generation of the `IMPACT-Z <https://github.com/impact-lbl/IMPACT-Z>`__ code.
ImpactX runs on modern GPUs or CPUs alike, provides user-friendly interfaces suitable for AI/ML workflows, has many :ref:`benchmarks <usage-examples>` to ensure its correctness, and an extensive documentation.

As a beam dynamics code, ImpactX uses the reference trajectory path length :math:`s` as the independent variable of motion to achieve large speedups.
The code includes the effects of externally applied fields from magnets and accelerating cavities as well as the effect of self-fields (space charge fields, CSR, wakefields, ...).
All particle tracking models are symplectic, and space charge is included by solving the Poisson equation in the beam rest frame.
The code may be used to model the dynamics of beams in both linear and ring accelerators.
See our :ref:`theory chapter <theory-concepts>` for details on our models, assumptions and concepts.

ImpactX is part of the `Beam, Plasma & Accelerator Simulation Toolkit (BLAST) <https://blast.lbl.gov>`__.
Please see this page for more information to select the code most appropriate for your application.


.. _contact:

Contact us
^^^^^^^^^^

We organize support as an open community, helping each other.

The `ImpactX GitHub repository <https://github.com/BLAST-ImpactX/impactx>`__ is our open communication platform.
Have a look at the action icons on the top right of the web page: feel free to watch the repo if you want to receive updates, or to star the repo to support the project.
For bug reports or to request new features, you can also open a new `issue <https://github.com/BLAST-ImpactX/impactx/issues>`__.

If you are starting to use ImpactX and have questions not answered by this documentation, have a look on our `discussion page <https://github.com/BLAST-ImpactX/impactx/discussions>`__.
There, you can find already answered questions, add your own, get help with installation procedures, help others and share ideas.

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   section#installation,
   section#usage,
   section#data-analysis,
   section#theory,
   section#development,
   section#maintenance,
   section#epilogue {
       display:none;
   }
   </style>

.. toctree::
   :hidden:

   coc
   acknowledge_us

Installation
------------
.. toctree::
   :caption: INSTALLATION
   :maxdepth: 1
   :hidden:

   install/users
   install/cmake
   install/hpc
   install/upgrade

Usage
-----
.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :hidden:

   usage/how_to_run
   usage/examples
   usage/python
   usage/parameters
   usage/dashboard
   usage/howto

Data Analysis
-------------
.. toctree::
   :caption: DATA ANALYSIS
   :maxdepth: 1
   :hidden:

   dataanalysis/dataanalysis

Theory
------
.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/concepts
   theory/coordinates_units
   theory/bending_elements
   theory/softedge_elements
   theory/beam_distribution
   theory/assumptions

Development
-----------
.. toctree::
   :caption: DEVELOPMENT
   :maxdepth: 1
   :hidden:

   developers/contributing
   developers/testing
   developers/documentation
   developers/repo_organization
   developers/implementation
   developers/doxygen
   developers/python
   developers/debugging
   developers/how_to_develop_with_llms

Maintenance
-----------
.. toctree::
   :caption: MAINTENANCE
   :maxdepth: 1
   :hidden:

   maintenance/release
..   maintenance/performance_tests

Epilogue
--------
.. toctree::
   :caption: EPILOGUE
   :maxdepth: 1
   :hidden:

   glossary
   acknowledgements
