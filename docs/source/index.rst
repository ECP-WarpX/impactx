:orphan:

ImpactX
-------

ImpactX is an s-based beam dynamics code including space charge effects.
This is the next generation of the `IMPACT-Z <https://github.com/impact-lbl/IMPACT-Z>`__ code.

Please contact us with any questions on it or if you like to contribute to its development.

.. _contact:

Contact us
^^^^^^^^^^

If you are starting using ImpactX, or if you have a user question, please pop in our `discussions page <https://github.com/ECP-WarpX/impactx/discussions>`__ and get in touch with the community.

The `ImpactX GitHub repo <https://github.com/ECP-WarpX/impactx>`__ is the main communication platform.
Have a look at the action icons on the top right of the web page: feel free to watch the repo if you want to receive updates, or to star the repo to support the project.
For bug reports or to request new features, you can also open a new `issue <https://github.com/ECP-WarpX/impactx/issues>`__.

We also have a `discussion page <https://github.com/ECP-WarpX/impactx/discussions>`__ on which you can find already answered questions, add new questions, get help with installation procedures, discuss ideas or share comments.

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   div#installation.section,
   div#usage.section,
   div#theory.section,
   div#data-analysis.section,
   div#development.section,
   div#maintenance.section,
   div#epilogue.section {
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
..   install/changelog
..   install/upgrade

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
   usage/workflows

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
