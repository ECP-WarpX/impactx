from __future__ import annotations

import amrex.space3d.amrex_3d_pybind
import impactx.impactx_pybind

__all__ = [
    "alpha",
    "convolve_fft",
    "deposit_charge",
    "derivative_charge",
    "unit_step",
    "w_l_csr",
    "w_l_rf",
    "w_t_rf",
]

def alpha(arg0: float) -> float:
    """
    Alpha Function
    """

def convolve_fft(
    arg0: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg1: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg2: float,
) -> amrex.space3d.amrex_3d_pybind.PODVector_real_std:
    """
    FFT Convolution
    """

def deposit_charge(
    arg0: impactx.impactx_pybind.ImpactXParticleContainer,
    arg1: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg2: float,
    arg3: float,
    arg4: bool,
) -> None:
    """
    Deposit Charge Distribution Function
    """

def derivative_charge(
    arg0: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg1: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg2: float,
    arg3: bool,
) -> None:
    """
    Derivative of Charge Profile Function
    """

def unit_step(arg0: float) -> float:
    """
    Step Function
    """

def w_l_csr(arg0: float, arg1: float, arg2: float) -> float:
    """
    CSR Wake Function
    """

def w_l_rf(arg0: float, arg1: float, arg2: float, arg3: float) -> float:
    """
    Longitudinal Resistive Wall Wake Function
    """

def w_t_rf(arg0: float, arg1: float, arg2: float, arg3: float) -> float:
    """
    Transverse Resistive Wall Wake Function
    """
