#!/bin/bash -l

# Copyright 2021-2023 ImpactX
#
# This file is part of ImpactX.
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 2
#SBATCH -J ImpactX
#SBATCH -A <proj>
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH -o ImpactX.o%j
#SBATCH -e ImpactX.e%j

# executable & inputs file or python interpreter & PICMI script here
EXE=./impactx
INPUTS=inputs_small

# each CPU node on Perlmutter (NERSC) has 64 hardware cores with
# 2x Hyperthreading/SMP
# https://en.wikichip.org/wiki/amd/epyc/7763
# https://www.amd.com/en/products/cpu/amd-epyc-7763
# Each CPU is made up of 8 chiplets, each sharing 32MB L3 cache.
# This will be our MPI rank assignment (2x8 is 16 ranks/node).

# threads for OpenMP and threaded compressors per MPI rank
export SRUN_CPUS_PER_TASK=16  # 8 cores per chiplet, 2x SMP
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun --cpu-bind=cores \
  ${EXE} ${INPUTS} \
  > output.txt
