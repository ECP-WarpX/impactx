from . import impactx_pybind as cxx
from .impactx_pybind import *  # noqa
from .madx_to_impactx import read_beam, read_lattice  # noqa

__version__ = cxx.__version__
__doc__ = cxx.__doc__
__license__ = cxx.__license__
__author__ = cxx.__author__

# at this place we can enhance Python classes with additional methods written
# in pure Python or add some other Python logic

# MAD-X reader for beamline lattice elements
#   adds an overload to existing methods
elements.KnownElementsList.load_file = lambda madx_file, nslice=1 : read_lattice(madx_file, nslice)  # noqa
