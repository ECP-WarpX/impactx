import os

# Python 3.8+ on Windows: DLL search paths for dependent
# shared libraries
# Refs.:
# - https://github.com/python/cpython/issues/80266
# - https://docs.python.org/3.8/library/os.html#os.add_dll_directory
if os.name == "nt":
    # add anything in the current directory
    pwd = __file__.rsplit(os.sep, 1)[0] + os.sep
    os.add_dll_directory(pwd)
    # add anything in PATH
    paths = os.environ.get("PATH", "")
    for p in paths.split(";"):
        p_abs = os.path.abspath(os.path.expanduser(os.path.expandvars(p)))
        if os.path.exists(p_abs):
            os.add_dll_directory(p_abs)

# import core bindings to C++
from . import impactx_pybind as cxx
from .extensions.ImpactXParIter import register_ImpactXParIter_extension
from .extensions.ImpactXParticleContainer import (
    register_ImpactXParticleContainer_extension,
)
from .impactx_pybind import *  # noqa
from .madx_to_impactx import read_beam, read_lattice  # noqa

__version__ = cxx.__version__
__doc__ = cxx.__doc__
__license__ = cxx.__license__
__author__ = cxx.__author__

# at this place we can enhance Python classes with additional methods written
# in pure Python or add some other Python logic

# MAD-X file reader for beamline lattice elements
elements.KnownElementsList.load_file = lambda self, madx_file, nslice=1: self.extend(
    read_lattice(madx_file, nslice)
)  # noqa

# MAD-X file reader for reference particle
RefPart.load_file = read_beam  # noqa

# Pure Python extensions to ImpactX types
register_ImpactXParIter_extension(cxx)
register_ImpactXParticleContainer_extension(cxx.ImpactXParticleContainer)

# include dashboard
from .dashboard import *  # noqa