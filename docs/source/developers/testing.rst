.. _developers-testing:

Testing
=======

Preparation
-----------

Prepare for running tests of ImpactX by :ref:`building ImpactX from source <install-developers>`.

In order to run our tests, you need to have a few :ref:`Python packages installed <install-dependencies>`:

.. code-block:: sh

   python3 -m pip install -U pip setuptools wheel pytest
   python3 -m pip install -r examples/requirements.txt

Run
---

You can run all our tests with:

.. code-block:: sh

   ctest --test-dir build --output-on-failure

Further Options
---------------

* help: ``ctest --test-dir build --help``
* list all tests: ``ctest --test-dir build -N``
* only run tests that have "FODO" in their name: ``ctest --test-dir build -R FODO``
