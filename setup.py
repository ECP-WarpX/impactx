#!/usr/bin/env python3
#
# Copyright 2021-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
import os
import platform
import re
import shutil
import subprocess
import sys

from setuptools import Extension, setup
from setuptools.command.build import build
from setuptools.command.build_ext import build_ext


class CopyPreBuild(build):
    def initialize_options(self):
        build.initialize_options(self)
        # We just overwrite this because the default "build" (and "build/lib")
        # clashes with directories many developers have in their source trees;
        # this can create confusing results with "pip install .", which clones
        # the whole source tree by default
        self.build_base = os.path.join("_tmppythonbuild", "impactx")

    def run(self):
        # remove existing build directory
        #   by default, this stays around. we want to make sure generated
        #   files like libimpactx.(so|pyd) are always only the
        #   ones we want to package and not ones from an earlier wheel's stage
        if os.path.exists(self.build_base):
            shutil.rmtree(self.build_base)

        # call superclass
        build.run(self)

        # copy Python module artifacts and sources
        dst_path = os.path.join(self.build_lib, "impactx")
        shutil.copytree(PYIMPACTX_libdir, dst_path, dirs_exist_ok=True)


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        from packaging.version import parse

        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake 3.24.0+ must be installed to build the following "
                + "extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        cmake_version = parse(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
        if cmake_version < parse("3.24.0"):
            raise RuntimeError("CMake >= 3.24.0 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + os.path.join(extdir, "impactx"),
            "-DCMAKE_VERBOSE_MAKEFILE=ON",
            "-DCMAKE_PYTHON_OUTPUT_DIRECTORY=" + extdir,
            "-DPython_EXECUTABLE=" + sys.executable,
            ## variants
            "-DImpactX_COMPUTE=" + ImpactX_COMPUTE,
            "-DImpactX_FFT:BOOL=" + ImpactX_FFT,
            "-DImpactX_MPI:BOOL=" + ImpactX_MPI,
            "-DImpactX_PRECISION=" + ImpactX_PRECISION,
            #'-DImpactX_PARTICLES_PRECISION=' + ImpactX_PARTICLES_PRECISION,
            "-DImpactX_PYTHON:BOOL=ON",
            ## dependency control (developers & package managers)
            #'-DImpactX_pyamrex_internal=' + ImpactX_pyamrex_internal,
            #'-DImpactX_pyamrex_repo=' + ImpactX_pyamrex_repo,
            #'-DImpactX_pyamrex_branch=' + ImpactX_pyamrex_branch,
            # PEP-440 conformant version from package
            "-DpyImpactX_VERSION_INFO=" + self.distribution.get_version(),
            #        see PICSAR and openPMD below
            ## static/shared libs
            "-DBUILD_SHARED_LIBS:BOOL=" + BUILD_SHARED_LIBS,
            ## Unix: rpath to current dir when packaged
            "-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON",
            "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=OFF",
            # Windows: has no RPath concept, all `.dll`s must be in %PATH%
            #          or same dir as calling executable
        ]
        # further dependency control (developers & package managers)
        if ImpactX_amrex_src:
            cmake_args.append("-DImpactX_amrex_src=" + ImpactX_amrex_src)
        # if ImpactX_pyamrex_src:
        #    cmake_args.append('-DImpactX_pyamrex_src=' + ImpactX_pyamrex_src)

        if CMAKE_INTERPROCEDURAL_OPTIMIZATION is not None:
            cmake_args.append(
                "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION="
                + CMAKE_INTERPROCEDURAL_OPTIMIZATION
            )
        if sys.platform == "darwin":
            cmake_args.append("-DCMAKE_INSTALL_RPATH=@loader_path")
        else:
            # values: linux*, aix, freebsd, ...
            #   just as well win32 & cygwin (although Windows has no RPaths)
            cmake_args.append("-DCMAKE_INSTALL_RPATH=$ORIGIN")

        cmake_args += extra_cmake_args

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(
                    cfg.upper(), os.path.join(extdir, "impactx")
                ),
            ]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

        # this environment variable is standardized in CMake
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # optional -j parameter in the build_ext call, not supported by pip
            if hasattr(self, "parallel") and self.parallel:
                build_args += ["-j{}".format(self.parallel)]
            else:
                build_args += ["-j2"]

        build_dir = os.path.join(self.build_temp, cfg)
        os.makedirs(build_dir, exist_ok=True)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_dir)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_dir)
        # note that this does not call install;
        # we pick up artifacts directly from the build output dirs


with open("./README.md", encoding="utf-8") as f:
    long_description = f.read()

# Allow to control options via environment vars.
#   Work-around for https://github.com/pypa/setuptools/issues/1712
# Pick up existing ImpactX libraries or...
PYIMPACTX_libdir = os.environ.get("PYIMPACTX_LIBDIR")

# ... build ImpactX libraries with CMake
#   note: changed default for SHARED, MPI, TESTING and EXAMPLES
ImpactX_COMPUTE = os.environ.get("IMPACTX_COMPUTE", "OMP")
ImpactX_FFT = os.environ.get("IMPACTX_FFT", "OFF")
ImpactX_MPI = os.environ.get("IMPACTX_MPI", "OFF")
ImpactX_PRECISION = os.environ.get("IMPACTX_PRECISION", "DOUBLE")
#   already prepared as a list 1;2;3
ImpactX_SPACEDIM = os.environ.get("IMPACTX_SPACEDIM", "3")
BUILD_SHARED_LIBS = os.environ.get("IMPACTX_BUILD_SHARED_LIBS", "OFF")
CMAKE_INTERPROCEDURAL_OPTIMIZATION = os.environ.get(
    "CMAKE_INTERPROCEDURAL_OPTIMIZATION", None
)
# CMake dependency control (developers & package managers)
ImpactX_amrex_src = os.environ.get("IMPACTX_AMREX_SRC")
ImpactX_amrex_internal = os.environ.get("IMPACTX_AMREX_INTERNAL", "ON")
ImpactX_amrex_repo = os.environ.get(
    "IMPACTX_AMREX_REPO", "https://github.com/AMReX-Codes/amrex.git"
)
ImpactX_amrex_branch = os.environ.get("IMPACTX_AMREX_BRANCH")
# ImpactX_pyamrex_src = os.environ.get('IMPACTX_PYAMREX_SRC')
# ImpactX_pyamrex_internal = os.environ.get('IMPACTX_PYAMREX_INTERNAL', 'ON')
# ImpactX_pyamrex_repo = os.environ.get('IMPACTX_PYAMREX_REPO',
#    'https://github.com/AMReX-Codes/pyamrex.git')
# ImpactX_pyamrex_branch = os.environ.get('IMPACTX_PYAMREX_BRANCH')

# extra CMake arguments
extra_cmake_args = []
for k, v in os.environ.items():
    extra_cmake_args_prefix = "IMPACTX_CMAKE_"
    if k.startswith(extra_cmake_args_prefix) and len(k) > len(extra_cmake_args_prefix):
        extra_cmake_args.append(
            "-D{0}={1}".format(k[len(extra_cmake_args_prefix) :], v)
        )

# https://cmake.org/cmake/help/v3.0/command/if.html
if ImpactX_MPI.upper() in ["1", "ON", "TRUE", "YES"]:
    ImpactX_MPI = "ON"
else:
    ImpactX_MPI = "OFF"

# for CMake
cxx_modules = []  # values: impactx or empty
cmdclass = {}  # build extensions

# externally pre-built: pick up pre-built pyImpactX libraries
if PYIMPACTX_libdir:
    cmdclass = dict(build=CopyPreBuild)
# CMake: build pyImpactX ourselves
else:
    cmdclass = dict(build_ext=CMakeBuild)
    cxx_modules.append(CMakeExtension("impactx"))

# Get the package requirements from the requirements.txt file
install_requires = []
with open("./requirements.txt") as f:
    install_requires = [line.strip("\n") for line in f.readlines()]
    if ImpactX_MPI == "ON":
        install_requires.append("mpi4py>=2.1.0")

# keyword reference:
#   https://packaging.python.org/guides/distributing-packages-using-setuptools
setup(
    name="impactx",
    # note PEP-440 syntax: x.y.zaN but x.y.z.devN
    version="24.08",
    packages=["impactx"],
    # Python sources:
    package_dir={"": "src/python"},
    # ImpactX authors:
    author="Axel Huebl, Chad Mitchell, Ryan Sandberg, Marco Garten, Ji Qiang, et al.",
    author_email="axelhuebl@lbl.gov, chadmitchell@lbl.gov, rsandberg@lbl.gov, mgarten@lbl.gov, jqiang@lbl.gov",
    # wheel/pypi packages:
    maintainer="Axel Huebl, Chad Mitchell, Ji Qiang",
    maintainer_email="axelhuebl@lbl.gov, chadmitchell@lbl.gov, jqiang@lbl.gov",
    description="ImpactX: the next generation of the IMPACT-Z beam dynamics code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords=(
        "research simulation particle-in-cell gpu accelerator physics pic particle beam-dynamics"
    ),
    url="https://impactx.readthedocs.io",
    project_urls={
        "Documentation": "https://impactx.readthedocs.io",
        # 'Tutorials': 'https://impactx-codes.github.io/impactx/tutorials_html/',
        "Doxygen": "https://impactx.readthedocs.io/en/latest/_static/doxyhtml",
        "Source": "https://github.com/ECP-WarpX/impactx",
        "DOI (source)": "https://doi.org/10.5281/zenodo.6954922",
        "DOI (paper)": "https://doi.org/10.48550/arXiv.2208.02382",
        "Tracker": "https://github.com/ECP-WarpX/impactx/issues",
    },
    # CMake: self-built as extension module
    ext_modules=cxx_modules,
    cmdclass=cmdclass,
    zip_safe=False,
    python_requires=">=3.8",
    tests_require=["numpy", "pandas", "pytest", "scipy"],
    install_requires=install_requires,
    # cmdclass={'test': PyTest},
    # platforms='any',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Natural Language :: English",
        "Environment :: Console",
        "Environment :: GPU",
        "Environment :: GPU :: NVIDIA CUDA",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        (
            "License :: OSI Approved :: " "BSD License"
        ),  # TODO: use real SPDX: BSD-3-Clause-LBNL
    ],
    # new PEP 639 format
    license="BSD-3-Clause-LBNL",
    license_files=["LICENSE"],
    entry_points={
        "console_scripts": [
            "impactx-dashboard=impactx.dashboard.__main__:main",
        ],
    },
)
