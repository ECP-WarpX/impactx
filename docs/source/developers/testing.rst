.. _developers-testing:

Testing
=======

Preparation
-----------

Prepare for running tests of ImpactX by :ref:`building ImpactX from source <install-developers>`.

In order to run our tests, you need to have a few :ref:`Python packages installed <install-dependencies>`:

.. code-block:: sh

   python3 -m pip install -U pip
   python3 -m pip install -U build packaging setuptools[core] wheel pytest pytest-benchmark
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

Single Precision
----------------

CI also runs the tests in a single-precision build, which is an in-development effort and will at a later point be fully supported.
To reproduce it locally:

.. code-block:: sh

   CXXFLAGS="-march=x86-64-v3" cmake -S . -B build_sp \
       -DImpactX_FASTMATH=ON      \
       -DImpactX_FFT=ON           \
       -DImpactX_PRECISION=SINGLE \
       -DImpactX_PYTHON=ON        \
       -DImpactX_SIMD=ON
   cmake --build build_sp -j 4
   ctest --test-dir build_sp --output-on-failure

A check that is tighter than single-precision round-off branches on the precision of the build.
Python scripts and pytest tests check ``impactx.Config.precision``.
Analysis scripts that only read openPMD data check the read back floating point type, e.g., ``is_double = initial["position_x"].dtype == np.float64``.

Known single-precision accuracy problems are set to a tolerance large enough to pass and have a ``FIXME`` comment that states the observed deviation and links the issue that tracks it.
