#!/usr/bin/env bash
#
set -eu -o pipefail

# adjust this to your project
export proj=m3906_g

# required dependencies
module load cmake/3.22.0 # cmake/git-20210830
module load PrgEnv-nvidia PrgEnv-gnu
module load cudatoolkit
module load cray-python/3.9.7.1

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 --ntasks-per-node=4 -t 1:00:00 -q interactive -C gpu --gpu-bind=single:1 -c 32 -G 4 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=4 -t 0:30:00 -q interactive -C gpu --gpu-bind=single:1 -c 32 -G 4 -A $proj"

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# necessary to use CUDA-Aware MPI and run a job
export CRAY_ACCEL_TARGET=nvidia80

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# compiler environment hints
export CC=cc #$(which gcc)
export CXX=CC #$(which g++)
export FC=ftn # $(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}

export CFLAGS="${CFLAGS:-} -O3 -ffast-math"
export CXXFLAGS="${CXXFLAGS:-} -O3 -ffast-math"
export FCLAGS="${FCFLAGS:-} -O3 -ffast-math"

# Get source upfront: ImpactX version 22.08
# Repo: https://github.com/ECP-WarpX/impactx
# Prepare via:
#   git clone --branch 22.08 https://github.com/ECP-WarpX/impactx ~/src/impactx
cd ~/src/impactx

# Benchmark inputs
BENCH_PARAMS="beam.npart=100000000 algo.space_charge=false diag.enable=false"

# CUDA benchmarks
rm -rf build_pm_cuda
cmake -S . -B build_pm_cuda -DImpactX_COMPUTE=CUDA
cmake --build build_pm_cuda -j 16

cd build_pm_cuda/bin
srun -N 1 --ntasks-per-node=1 --threads-per-core=1 -t 0:30:00 -q interactive -C gpu \
    -c 64 -G 1 \
    -A $proj \
    ./impactx ../../examples/iota_lattice/input_iotalattice.in ${BENCH_PARAMS} \
    > output_a100_1.txt
srun -N 1 --ntasks-per-node=2 --threads-per-core=1 -t 0:30:00 -q interactive -C gpu \
    -c 32 -G 2 \
    -A $proj \
    ./impactx ../../examples/iota_lattice/input_iotalattice.in ${BENCH_PARAMS} \
    > output_a100_2.txt
cd -

# CPU benchmarks
export MPICH_GPU_SUPPORT_ENABLED=0
#export SLURM_CPU_BIND="cores"

rm -rf build_pm_noacc
cmake -S . -B build_pm_noacc -DImpactX_COMPUTE=NOACC
cmake --build build_pm_noacc -j 16

cd build_pm_noacc/bin
srun -N 1 --ntasks-per-node=1 --threads-per-core=1 -t 0:30:00 -q interactive -C gpu -c 2 \
    -A $proj \
    ./impactx ../../examples/iota_lattice/input_iotalattice.in ${BENCH_PARAMS} \
    > output_cpu_1.txt
srun -N 1 --ntasks-per-node=8 --threads-per-core=1 -t 0:30:00 -q interactive -C gpu -c 2 \
    -A $proj \
    ./impactx ../../examples/iota_lattice/input_iotalattice.in ${BENCH_PARAMS} \
    > output_cpu_8.txt
srun -N 1 --ntasks-per-node=16 --threads-per-core=1 -t 0:30:00 -q interactive -C gpu -c 2 \
    -A $proj \
    ./impactx ../../examples/iota_lattice/input_iotalattice.in ${BENCH_PARAMS} \
    > output_cpu_16.txt
srun -N 1 --ntasks-per-node=32 --threads-per-core=1 -t 0:30:00 -q interactive -C gpu -c 2 \
    -A $proj \
    ./impactx ../../examples/iota_lattice/input_iotalattice.in ${BENCH_PARAMS} \
    > output_cpu_32.txt
cd -
