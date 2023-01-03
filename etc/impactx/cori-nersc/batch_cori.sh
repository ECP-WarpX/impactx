#!/bin/bash -l

# Copyright 2019-2023 Maxence Thevenet
#
# This file is part of ImpactX.
#
# License: BSD-3-Clause-LBNL


#SBATCH -N 2
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -C knl
#SBATCH -S 4
#SBATCH -J <job name>
#SBATCH -A <allocation ID>
#SBATCH -e ImpactX.e%j
#SBATCH -o ImpactX.o%j

export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# KNLs have 4 hyperthreads max
export CORI_MAX_HYPETHREAD_LEVEL=4
# We use 64 cores out of the 68 available on Cori KNL,
# and leave 4 to the system (see "#SBATCH -S 4" above).
export CORI_NCORES_PER_NODE=64

# Typically use 8 MPI ranks per node without hyperthreading,
# i.e., OMP_NUM_THREADS=8
export IMPACTX_NMPI_PER_NODE=8
export IMPACTX_HYPERTHREAD_LEVEL=1

# Compute OMP_NUM_THREADS and the thread count (-c option)
export CORI_NHYPERTHREADS_MAX=$(( ${CORI_MAX_HYPETHREAD_LEVEL} * ${CORI_NCORES_PER_NODE} ))
export IMPACTX_NTHREADS_PER_NODE=$(( ${IMPACTX_HYPERTHREAD_LEVEL} * ${CORI_NCORES_PER_NODE} ))
export OMP_NUM_THREADS=$(( ${IMPACTX_NTHREADS_PER_NODE} / ${IMPACTX_NMPI_PER_NODE} ))
export IMPACTX_THREAD_COUNT=$(( ${CORI_NHYPERTHREADS_MAX} / ${IMPACTX_NMPI_PER_NODE} ))

# for async_io support: (optional)
export MPICH_MAX_THREAD_SAFETY=multiple

srun --cpu_bind=cores -n $(( ${SLURM_JOB_NUM_NODES} * ${IMPACTX_NMPI_PER_NODE} )) -c ${IMPACTX_THREAD_COUNT} \
  <path/to/executable> <input file> \
  > output.txt
