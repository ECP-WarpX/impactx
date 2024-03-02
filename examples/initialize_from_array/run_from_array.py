import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.constants import c, e, m_e, epsilon_0
import scipy.optimize as opt

try:
    import cupy as cp

    cupy_available = True
except ImportError:
    cupy_available = False

from openpmd_viewer.addons import LpaDiagnostics

from impactx import Config, ImpactX, RefPart, distribution, elements, ImpactXParIter
from impactx import coordinate_transformation, CoordSystem

import torch