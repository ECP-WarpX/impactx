"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Chad Mitchell, Axel Huebl
License: BSD-3-Clause-LBNL
"""


def register_RFCavity_extension(cls):
    """Extend RFCavity with an alternative constructor that accepts
    raw on-axis field data and computes Fourier coefficients internally.

    Parameters
    ----------
    cls : type
        The pybind11 ``RFCavity`` class to extend.
    """
    _original_init = cls.__init__

    def _new_init(
        self,
        *,
        ds,
        escale,
        freq,
        phase,
        cos_coefficients=None,
        sin_coefficients=None,
        z=None,
        field_on_axis=None,
        ncoef=None,
        z_data=None,
        wake_x_data=None,
        wake_y_data=None,
        wake_z_data=None,
        dx=0,
        dy=0,
        rotation=0,
        aperture_x=0,
        aperture_y=0,
        mapsteps=10,
        nslice=1,
        name=None,
    ):
        has_coefficients = cos_coefficients is not None or sin_coefficients is not None
        has_field_data = z is not None or field_on_axis is not None or ncoef is not None

        if has_coefficients and has_field_data:
            raise ValueError(
                "RFCavity: provide either (cos_coefficients, sin_coefficients) "
                "or (z, field_on_axis, ncoef), not both."
            )

        if has_field_data:
            if z is None or field_on_axis is None or ncoef is None:
                raise ValueError(
                    "RFCavity: when using field data, all three parameters "
                    "'z', 'field_on_axis', and 'ncoef' must be provided."
                )
            from ..fourier import fourier_coefficients

            cos_coefficients, sin_coefficients = fourier_coefficients(
                z, field_on_axis, ncoef
            )

        kwargs = dict(
            ds=ds,
            escale=escale,
            freq=freq,
            phase=phase,
            dx=dx,
            dy=dy,
            rotation=rotation,
            aperture_x=aperture_x,
            aperture_y=aperture_y,
            mapsteps=mapsteps,
            nslice=nslice,
        )
        if name is not None:
            kwargs["name"] = name
        if cos_coefficients is not None:
            kwargs["cos_coefficients"] = cos_coefficients
        if sin_coefficients is not None:
            kwargs["sin_coefficients"] = sin_coefficients
        if z_data is not None:
            kwargs["z_data"] = z_data
        if wake_x_data is not None:
            kwargs["wake_x_data"] = wake_x_data
        if wake_y_data is not None:
            kwargs["wake_y_data"] = wake_y_data
        if wake_z_data is not None:
            kwargs["wake_z_data"] = wake_z_data

        _original_init(self, **kwargs)

    _new_init.__doc__ = _original_init.__doc__
    cls.__init__ = _new_init
