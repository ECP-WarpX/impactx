"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def ix_pc_to_df(self, local=True, comm=None, root_rank=0):
    """
    Copy all particles into a pandas.DataFrame

    Parameters
    ----------
    self : ImpactXParticleContainer_*
        The particle container class in ImpactX
    local : bool
        MPI-local particles
    comm : MPI Communicator
        if local is False, this defaults to mpi4py.MPI.COMM_WORLD
    root_rank : MPI root rank to gather to
        if local is False, this defaults to 0

    Returns
    -------
    A concatenated pandas.DataFrame with particles from all levels.

    Returns None if no particles were found.
    If local=False, then all ranks but the root_rank will return None.
    """
    df = super(type(self), self).to_df(local=local, comm=comm, root_rank=root_rank)

    # rename columns according to our attribute names
    if df is not None:
        # todo: check if currently in fixed s or fixed t and pick name accordingly

        names = []
        for n in self.RealAoS_names:
            names.append(n)
        names.append("cpuid")
        for n in self.RealSoA_names:
            names.append(n)

        df.columns.values[0 : len(names)] = names

    return df


def register_ImpactXParticleContainer_extension(ixpc):
    """ImpactXParticleContainer helper methods"""

    # register member functions for ImpactXParticleContainer
    ixpc.to_df = ix_pc_to_df
