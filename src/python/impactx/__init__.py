from . import impactx_pybind
from .MADXParser import *  # noqa
from .impactx_pybind import *  # noqa
from .madx_to_impactx import madx2impactx_beam, madx2impactx_lattice  # noqa

__version__ = impactx_pybind.__version__
__doc__ = impactx_pybind.__doc__
__license__ = impactx_pybind.__license__
__author__ = impactx_pybind.__author__

# at this place we can enhance Python classes with additional methods written
# in pure Python or add some other Python logic
#
