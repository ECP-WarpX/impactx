.. _usage-python:

Dashboard
==================

This document provides instructions on how to use the `impactx-dashboard`.

Installing Dependencies
-----------------------
From the `impactx/src/python/impactx/dashboard` directory, execute:

.. code-block:: bash

   pip install -r requirements.txt

Launching the Dashboard
-----------------------

There are two ways to launch the Dashboard GUI, depending on your use case:

- **Using the CLI Command**: This is the recommended method for most users. It is straightforward and requires minimal setup.
- **Using the Python Module**: This method is more flexible and is particularly useful for developers who need to make direct
changes to the dashboard, integrate it with other Python scripts, or run it in a specific environment.

1. **Direct Command Execution**
   From the base ImpactX directory, execute:

   .. code-block:: bash

      impactx-dashboard

2. **Python Module Execution**
   From the `impactx/src/python/impactx` directory, execute:

   .. code-block:: bash

      python -m dashboard

3. **Access the GUI**

   Open a web browser and navigate to `http://localhost:8080` to access the Dashboard GUI.

   #TODO accessing through Jupyter
