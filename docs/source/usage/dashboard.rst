.. _usage-dashboard:

Dashboard
=========

ImpactX Dashboard is a browser-based interface to ImpactX.
It provides a graphical interface to a subset of ImpactX functionality.

.. note::

   ImpactX Dashboard is provided as a preview and continues to be developed.
   Thus, it may not yet include all the features available in ImpactX.
   Let us know in GitHub `discussions <https://github.com/ECP-WarpX/impactx/discussions>`__ and `issues <https://github.com/ECP-WarpX/impactx/issues>`__ how it works for you and what kind of workflows you would like to run in it.


Launching the Dashboard
-----------------------

The ImpactX Dashboard can be run in two modes, as a standalone browser application or inside a Jupyter notebook.

1. **Standalone browser application:**
   After installation of ImpactX including the Python modules, launch:

   .. code-block:: bash

      impactx-dashboard

2. **JupyterLab:**
   Start `JupyterLab <https://jupyter.org/install>`__, e.g., logging into a Jupyter service at an HPC center or locally on your computer using:

   .. code-block:: bash

      jupyter-lab

   Inside JupyterLab, run the following Python code in a notebook to initialize and display the dashboard:

   .. code-block:: python

      from impactx.dashboard import JupyterApp

      # Create new application instance
      app = JupyterApp()

      # Start the server and wait for the UI to be ready
      await app.ui.ready

      # Display the UI in the JupyterLab notebook
      app.ui


Navigation
----------

The user-friendly interface includes multiple tabs and menu options, intended to be navigated from top to bottom:

- **Input Tab**: Allows to adjust simulation input parameters.
- **Run Tab**: Enables to run simulations and monitor their progress.
- **Analyze Tab**: Provides tools to visualize and analyze simulation results.

.. figure:: https://gist.githubusercontent.com/ax3l/b56aa3c3261f9612e276f3198b34f771/raw/11bfe461a24e1daa7fd2d663c686b0fcc2b6e305/dashboard.png
   :align: center
   :width: 75%
   :alt: phase space ellipse

   Input section in the dashboard.


Developers
----------

Additional Dependencies
"""""""""""""""""""""""

Additional dependencies to ImpactX for the dashboard are included relative ImpactX source directory:

.. code-block:: bash

   python -m pip install -r src/python/impactx/dashboard/requirements.txt

Python Module
"""""""""""""

After installing only the ImpactX Python bindings, one can directly run the dashboard modules from the source tree during development, too.
For this, navigate in the ImpactX source directory to the ``src/python/impactx`` directory and run:

.. code-block:: bash

   python -m dashboard
