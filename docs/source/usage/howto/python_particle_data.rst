.. _usage-howto-python-particle-data:

Accessing Beam Particle Data
============================

Selecting the beam
------------------

The beam is accessed through :py:attr:`sim.beam <impactx.ImpactX.beam>`, where ``sim`` is the
:py:class:`~impactx.ImpactX` instance set up as described in :ref:`usage-howto-python-extend-run`.
This returns an :py:class:`impactx.ParticleContainer`, distributed over MPI ranks.

.. code-block:: python

   pc = sim.beam

The available per-particle attributes are:

* **Positions** relative to the reference particle: ``position_x``, ``position_y``, ``position_t`` (all in meters; ``position_t`` is :math:`c \Delta t`).
* **Momenta** relative to and normalized by the reference momentum: ``momentum_x``, ``momentum_y``, ``momentum_t``.
* **Spin** components (only if ``sim.spin = True``): ``spin_x``, ``spin_y``, ``spin_z``.
* **Charge over mass** ``qm`` (in 1/eV), which is currently inconsistently used in ImpactX.
* **Macroparticle weight** ``w`` (unitless), how many physical particles a simulation particle represents.
* **Unique ID** ``idcpu``, a 64bit integer that is unique to a particle over the lifetime of a simulation.

The :py:attr:`~impactx.ParticleContainer.ref` attribute provides the reference particle (a
:py:class:`impactx.RefPart`). All phase-space coordinates above are relative to the reference particle.

.. note::

   The independent variable in ImpactX is the reference-trajectory path length ``s`` (not time).
   See :ref:`theory-coordinates-and-units` for the full coordinate and unit conventions.

Accessing/modifying the underlying particle data
------------------------------------------------

There are two ways to access particle data, with different trade-offs between convenience and performance.
(For an in-depth discussion of the underlying pyAMReX API, see the
`pyAMReX particles guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__.)

.. tab-set::

   .. tab-item:: Global access through a pandas DataFrame (read-only)

      The :py:meth:`~impactx.ParticleContainer.to_df` method returns a
      `pandas DataFrame <https://pandas.pydata.org/docs/user_guide/dsintro.html#dataframe>`__
      containing the particle data. The columns of the DataFrame are the particle attributes
      (e.g. ``position_x``, ``momentum_x``, ``w``, ``id``), and each row corresponds to one
      macroparticle across all tiles on the current MPI rank.

      .. warning::

         The DataFrame is a *copy* of the particle data. Modifying it does not write back
         to the simulation.

      .. note::

         ``to_df`` is convenient because it concatenates all particles across tiles
         into a single table. However, it incurs copies and, on GPU runs, CPU↔GPU data transfers.
         It is well-suited to debugging and visualization, and not to performance-critical paths.

      .. code-block:: python

         pc = sim.beam

         # local particles only (default): returns particles on the current MPI rank
         df = pc.to_df(local=True)
         if df is not None:
             print("Available attributes:", list(df.columns))
             print("Number of particles:", len(df))
             print("position_x:", df["position_x"])

         # positions in the lab frame, by adding the reference particle position
         ref = pc.ref
         x_lab = ref.x + df["position_x"]

      To gather particles from *all* MPI ranks onto the root rank, pass ``local=False``:

      .. code-block:: python

         df_global = pc.to_df(local=False)
         # df_global is non-None only on the root rank

   .. tab-item:: Explicit loop over particle tiles

      The :py:class:`impactx.ImpactXParIter` iterator gives direct, zero-copy access to the
      Struct-of-Arrays storage of each tile on each mesh-refinement level. This avoids the copies
      performed by ``to_df`` and is the right choice for performance-critical analysis or in-place
      modification of the beam.

      .. code-block:: python

         from impactx import ImpactXParIter

         pc = sim.beam

         for lvl in range(pc.finest_level + 1):
             for pti in ImpactXParIter(pc, level=lvl):
                 soa = pti.soa().to_xp()   # NumPy (CPU) or CuPy (GPU)
                 x  = soa.real["position_x"]
                 px = soa.real["momentum_x"]

                 # in-place modification works:
                 x[:] += 1.0e-6   # shift every particle by 1 µm

      Inside a :py:class:`impactx.elements.Programmable` ``beam_particles`` callback, the
      framework hands you a single ``pti`` directly — no need to construct the iterator yourself.
      See :ref:`usage-howto-python-extend-callbacks` for a full example.

Reference particle
------------------

The reference particle is accessed via :py:attr:`pc.ref <impactx.ParticleContainer.ref>` and
returns a :py:class:`impactx.RefPart`. It exposes the integrated path length ``s``, the lab-frame
position/momentum (``x``, ``y``, ``z``, ``px``, ``py``, ``pz``), the time-of-flight ``t``, the
energy-related quantities (``pt``, ``beta``, ``gamma``, ``beta_gamma``), and the species mass/charge.

.. code-block:: python

   ref = sim.beam.ref
   print(f"s = {ref.s:.3f} m, gamma = {ref.gamma:.3f}")

When using :py:class:`~impactx.elements.Programmable` with separate ``ref_particle`` and
``beam_particles`` callbacks, the reference particle is pushed *before* the beam particles,
so inside the ``beam_particles`` callback ``refpart`` already reflects the updated reference
state for the current slice.

Adding new particles
--------------------

New particles can be appended to the beam at any time after :py:meth:`~impactx.ImpactX.init_grids`
using :py:meth:`~impactx.ParticleContainer.add_n_particles`. See its API reference for parameters.

.. _usage-howto-python-particle-data-fixed-t:

Exchanging particles with time-based codes
------------------------------------------

Codes that use time as the independent variable (e.g., WarpX) describe the beam at fixed time :math:`t`, while ImpactX describes it at a fixed reference position :math:`s`.
To read, modify, or add particles at fixed :math:`t`, work inside a :py:meth:`~impactx.ParticleContainer.at_fixed_t` block:

.. code-block:: python

   with sim.beam.at_fixed_t() as beam:
       for lvl in range(beam.finest_level + 1):
           for pti in ImpactXParIter(beam, level=lvl):
               soa = pti.soa().to_xp()
               z = soa.real["position_t"]   # holds z at fixed t
               pz = soa.real["momentum_t"]  # holds pz at fixed t
               ...  # hand x, y, z, px, py, pz to the other code, write back its result

   # here, the beam is at fixed s again

Inside the block, the arrays keep their fixed-:math:`s` names, but ``position_t`` and ``momentum_t`` hold the longitudinal position :math:`z` (in meters) and momentum deviation :math:`p_z`, both relative to the reference particle.
When the block ends, normally or through an exception, the beam is transformed back to fixed :math:`s`.

Each transformation uses the reference energy ``beam.ref.pt`` at the moment it runs: the transformation to fixed :math:`t` uses the reference particle as the block starts, the transformation back to fixed :math:`s` uses the reference particle as the block ends.
If the other code changes the beam energy, update the reference particle inside the block (e.g., with :py:meth:`~impactx.RefPart.set_kin_energy_MeV`) and write the particle coordinates relative to the updated reference particle.

The :ref:`initialize from arrays <examples-from-array>` and :ref:`ML surrogate <examples-ml-surrogate>` examples use this pattern.

See also
--------

* :ref:`usage-howto-python-extend`: where and when to install callbacks.
* :ref:`usage-howto-python-field-data`: how to access the space-charge fields.
* :ref:`dataanalysis`: offline analysis from the openPMD output of the :py:class:`~impactx.elements.BeamMonitor`.
* `pyAMReX particles guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles>`__.
