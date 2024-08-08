.. _usage-dashboard:

Dashboard
=========

This document provides instructions on how to launch the ImpactX Dashboard, a browser-based interface to ImpactX.


Launching the Dashboard
-----------------------

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


Navigating the GUI
------------------

The GUI is designed with a user-friendly interface that includes multiple tabs and menu options:

- **Input Tab**: Allows to adjust simulation input parameters.
- **Run Tab**: Enables to run simulations and monitor their progress.
- **Analyze Tab**: Provides tools to visualize and analyze simulation results.

.. figure:: https://gist.githubusercontent.com/ax3l/b56aa3c3261f9612e276f3198b34f771/raw/25a9e2709a449fb4fd1c2ffb6c961716e70d8b32/dashboard.png
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
For this, navigate to the ``impactx/src/python/impactx/dashboard`` directory and run:

   .. code-block:: bash

      python -m dashboard
