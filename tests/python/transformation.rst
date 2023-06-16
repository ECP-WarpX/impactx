.. _tests-transformation:

Transformation
==============

Test the t/s transformations on an electron beam.

We use a long electron beam, :math:`L_z=1` cm, with significant correlations in :math:`x-px`, :math:`y-py`, and :math:`t-pt`.
The beam has average energy 1 GeV.

This tests that the t/s transforms are inverses of each other
Specifically, in this test the :math:`t`- and :math:`s`-coordinates of the beam must differ substantially
and the forward-inverse transformed coordinates must agree with the initial coordinates.
That is, we require that :math:`to_fixed_s` ( :math:`to_fixed_t` (initial beam)) = initial beam.


Run
---

This file is run from :ref:`pytest <developers-testing>`.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: test_transformation.py
          :language: python3
          :caption: You can copy this file from ``tests/python/test_transformation.py``.
