# AGENTS Guidelines for This Repository

## Project Overview

ImpactX is a high-performance beam dynamics code for particle accelerators with collective effects, the next generation of the IMPACT-Z code. It is built on top of AMReX (adaptive mesh refinement framework) and shares the ABLASTR abstraction layer with WarpX. Particle tracking uses the reference trajectory path length `s` as the independent variable, with symplectic integrators and self-fields (space charge, CSR, wakefields). It runs on CPU/OpenMP, CUDA, HIP, and SYCL backends.

## Domain Context

Reference `docs/source/glossary.rst` as context for WarpX terminology, acronyms, and physics concepts.

## Development Environment

If you cannot find the `cmake` or `ctest` command, activate the conda environment named `impactx-cpu-mpich-dev` before running shell commands that compile or test ImpactX.
The `conda` command may be called `conda`, `mamba`, or `micromamba` depending on the system (check aliases and available binaries):
```bash
conda activate impactx-cpu-mpich-dev
```
If this environment does not yet exist, create it as described in `docs/source/install/dependencies.rst`.

## Build Commands

The cmake build directory is always inside the repository root (or worktree root). Never look for or create a build directory outside of the current working directory.

```bash
# Configure (common development build with Python bindings)
cmake --fresh -S . -B build -DImpactX_PYTHON=ON

# Build
cmake --build build -j 8

# Install Python bindings
cmake --build build --target pip_install

# Faster rebuild (skip dependency checks)
cmake --build build -j 8 --target pip_install_nodeps

# Faster link times during development
cmake -S . -B build -DImpactX_PYTHON=ON -DImpactX_PYTHON_IPO=OFF -DpyAMReX_IPO=OFF
```

Key CMake options:
- `-DImpactX_COMPUTE=OMP` — backend: `NOACC`, `OMP`, `CUDA`, `SYCL`, `HIP`
- `-DImpactX_MPI=ON` — multi-node support
- `-DImpactX_PYTHON=ON` — Python bindings
- `-DImpactX_OPENPMD=ON` — HDF5/ADIOS I/O (on by default)
- `-DImpactX_FFT=ON` — FFT-based solvers (IGF space charge, CSR)
- `-DImpactX_PRECISION=DOUBLE` — floating point precision: `SINGLE` or `DOUBLE`
- `-DImpactX_SIMD=ON` — CPU SIMD acceleration

## Testing

Tests use CTest. Each example in `examples/` has an input file (or Python script), an analysis script, and optionally a plot script. Python unit tests live under `tests/python/` and are run via pytest.

You can run the Python unit tests directly with `pytest` from `tests/python/`, but this requires that `cmake --build build --target pip_install` has been run first so the Python package is installed and importable.
The pytest-based tests are also registered in CTest. Running them via `ctest` is often preferable during development because CTest sets the needed environment variables (especially `PYTHONPATH`) so the build tree is preferred automatically.

```bash
# List tests
ctest --test-dir build -N

# Run all tests
ctest --test-dir build --output-on-failure

# Run specific test by regex
ctest --test-dir build -R FODO

# Run exact test
ctest --test-dir build -R "^FODO\..*"
```

Test output goes to `build/bin/<test_name>/`.

- When debugging/fixing tests: do not modify the tolerance of assert statements in the Python analysis files just to make the tests pass (unless explicitly asked to do so).

### Adding a Test

Use `add_impactx_test()` in `examples/CMakeLists.txt`:

```cmake
add_impactx_test(<name>
    examples/<dir>/<input_or_script>   # input file (.in) or Python script (.py)
    OFF                                 # is_mpi: ON to run with mpirun
    examples/<dir>/analysis_<...>.py    # validation script
    examples/<dir>/plot_<...>.py        # optional plot script
)
```

Tests must be quick to run on a 2 core CI CPU runner (ideally <30sec) and be written in a CPU/GPU portable way.
Analysis functions validate the outputs are as expected (e.g., Twiss parameters, emittances, moments).

## Breaking Changes

A change that makes an existing user script behave differently, or stop working, is
documented in the upgrade guide `docs/source/install/upgrade.rst`, under a section for the
release it lands in (newest first).

Write what a user has to do:

- what changed, in terms of the API they call;
- what to look for in their own scripts;
- the replacement, as code.

This is the one place where "it used to do X" belongs. Everywhere else -- docstrings, the
rest of the manual -- documents current behavior only. There is no changelog in the docs;
per-release notes live in the GitHub release.

## Auto-Generated Files

- **`.pyi` stub files** under `src/python/impactx/`: Never modify. They are auto-generated on the `development` branch after a PR is merged (see `.github/workflows/stubs.yml`).

## Linting and Formatting

Pre-commit hooks handle formatting:
```bash
python -m pip install -U pre-commit
pre-commit install
pre-commit run -a  # run on all files
```

- **Python**: Ruff (configured in `pyproject.toml`)
- **C++**: no auto-formatter; follow the style below and existing code

Commits should limit any formatting changes of unchanged code.

## Code Architecture

### Source Layout

- `src/ImpactX.H/.cpp` — main simulation class `impactx::ImpactX`
- `src/main.cpp` — standalone executable entry point
- `src/tracking/` — tracking loops (`track_particles`, `track_envelope`, `track_reference`)
- `src/elements/` — lattice elements (drifts, dipoles, quadrupoles, RF cavities, plasma lenses, apertures, …); `All.H` lists the `KnownElements` variant
- `src/particles/` — `ImpactXParticleContainer`, distributions, integrators, charge deposition, collect-lost
- `src/envelope/` — linear envelope tracking via the covariance matrix
- `src/initialization/` — AMReX setup, input parsing, mesh refinement, distribution & element initialization
- `src/diagnostics/` — reduced beam characteristics, emittance invariants, openPMD output
- `src/python/` — pybind11 bindings (`pyImpactX.cpp`, element/distribution bindings)
- `src/python/impactx/` — Python module (dashboard, MADX import, plotting helpers)
- `examples/` — physics-driven integration tests (FODO, chicane, cyclotron, IOTA, …)
- `tests/python/` — pytest-based unit tests

### Key Patterns

- Core simulation state (AMR hierarchy, fields, particle containers) lives in `impactx::initialization::AmrCoreData`, accessed via `ImpactX::amr_data`
- Lattice stored as `std::list<elements::KnownElements> m_lattice` on `ImpactX`; new elements are added to the `KnownElements` variant in `src/elements/All.H`
- Space-charge and CSR solvers and many utilities come from `ABLASTR` (shared with WarpX, fetched into `build/_deps/fetchedablastr-src/`)
- Python bindings are generated via `pybind11` and `pyAMReX`; after C++ changes, rebuild with the `pip_install` target
- The independent variable of motion is the reference trajectory path length `s` (not time); all integrators are symplectic
- `amrex::ParallelFor` promises the compiler that loop iterations are independent (it applies a CPU SIMD pragma). Kernels where different iterations can write the same memory location (e.g. deposition, histogram bins, shared counters, etc.) must use `amrex::For` instead, and whole-loop sums/maxima the `amrex::Reduce` functions. `amrex::Gpu::Atomic` and `amrex::HostDevice::Atomic` operations are plain non-atomic updates on serial CPU builds and do not make a `ParallelFor` safe. See https://github.com/BLAST-WarpX/warpx/issues/7097

## C++ Style

- 4 spaces indentation, no tabs, max 100 chars/line
- Member variables prefixed with `m_`
- CamelCase for files and classes (e.g., `ImpactXParticleContainer.H`)
- Space before `()` in function declarations and definitions
- Opening brace `{` on new line for functions/classes
- Always use curly braces for single-statement blocks
- Use `amrex::` namespace prefix (avoid `using namespace amrex`); ImpactX code lives in `namespace impactx { … }`

### Include Order

In `.cpp` files: (1) corresponding header, (2) ImpactX headers, (3) ImpactX forward declarations, (4) ABLASTR headers, (5) AMReX headers, (6) AMReX forward declarations, (7) third-party headers (pybind11, openPMD, …), (8) standard library. Each group alphabetically sorted with blank lines between groups.

## Version Control

- Main branch: `development` (not `main`)
- Fork-and-branch workflow; PRs target `development`
- Pull requests with features and bug fixes need to add a test for coverage.
