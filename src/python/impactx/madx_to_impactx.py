from impactx import elements

def madx2impactx(parsed_beamline,nslice=25):
    """
    Function that converts a list of elements in the MADXParser former into ImpactX format
    :param parsed_beamline: list of dictionaries
    :return: list of translated dictionaries
    """

    madx_to_impactx_dict = {
        "MARKER": "None",
        "DRIFT": "Drift",
        "SBEND": "Sbend",  # Sector Bending Magnet
        "QUAD": "Quad",  # Quadrupole
        "DIPEDGE": "DipEdge",
        "MULTIPOLE": "Multipole",
        "NLLENS": "NonlinearLens",
        # TODO Figure out how to identify these
        "ShortRF": "ShortRF",
        "ConstF": "ConstF",
    }

    impactx_beamline = []

    for d in parsed_beamline:

        if d['type'] in [k.casefold() for k in list(madx_to_impactx_dict.keys())]:
            if d['name'] == "drift":
                impactx_beamline.append(
                    elements.Drift(ds=d['l'], nslice=nslice)
                )
            elif d['name'] == "quadrupole":
                impactx_beamline.append(
                    elements.Quad(ds=d['l'], k=d['k1'], nslice=nslice)
                )
        else:
            raise NotImplementedError(
                "The beamline element named ",
                d['name'],
                "of type ",
                d['type'],
                "is not implemented in impactx.elements.",
                "Available elements are:",
                list(madx_to_impactx_dict.keys())
            )
    return impactx_beamline


