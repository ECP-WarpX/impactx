#!/bin/bash -l

#SBATCH -A <project id>
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=impactx
#SBATCH --output=ImpactX.o%j
#SBATCH --error=ImpactX.e%j
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8

# executable & inputs file or ...
EXE=./impactx
INPUTS=inputs
# ... python interpreter & PICMI script here
#EXE=python3
#INPUTS=run_fodo.py

date

# note (12-12-22)
# this environment setting is currently needed on LUMI to work-around a
# known issue with Libfabric
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

# Seen since August 2023 seen on OLCF (not yet seen on LUMI?)
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors
#export MPICH_SMP_SINGLE_COPY_MODE=NONE
#export FI_CXI_RX_MATCH_MODE=software

# Seen since August 2023
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software

# The LUMI-C compute nodes are configured with 4 NUMA zones ("quadrant mode") with 2 CCDs per quadrant.
# https://docs.lumi-supercomputer.eu/hardware/lumic/
#
# LUMI documentation suggests using the following CPU bind
# see https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumic-job/
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun                 \
  ${EXE} ${INPUTS}   \
  | tee outputs.txt
