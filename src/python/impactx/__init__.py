from . import impactx_pybind as cxx
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
RefPart.load_file = lambda self, madx_file: read_beam(self, madx_file)  # noqa
