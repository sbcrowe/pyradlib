# -*- coding: utf-8 -*-
"""Xoft Axxent dose calculation module.

This module provides functions assisting in the determination of dose for the Xoft Axxent system.

Notes
-----
This module assumes that vendor supplied "Balloon Applicator Atlas.xml" and "xoft_axxent.xlsx" are 
provided in a "data" subdirectory.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = ["Rachael Wilks"]
__license__ = "GPL3"

# import required code
import math
import numpy as np
import os
import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.interpolate import interpn

# define global variables
_reference_air_kerma = 10000
_dose_rate_constant = 8.520

# find data location
_data_root_path = os.path.join(os.fspath(Path(__file__).parent.resolve()), "data")


def dwell_data(plan_name: str):
    """Imports template plan dwell positions from vendor-provided xml files.

    Parameters
    ----------
    plan_name : str
        Name of the template plan, for example '3-4_30cc'.

    Returns
    -------
    dwell_positions, dwell_times : array_like
        Lists of dwell positions and times for corresponding template plan.
    """
    balloon = plan_name[0:3]
    inflation = plan_name[4:]
    tree = ET.parse(os.path.join(_data_root_path, "Balloon Applicator Atlas.xml"))
    for parameter_1 in tree.iter("Parameter_1"):
        for atlas_parameter_1 in parameter_1.iter("AtlasParameter1"):
            if atlas_parameter_1.findall("Description")[0].text.startswith(
                balloon
            ):  # correct balloon
                for parameter_2 in atlas_parameter_1.iter("Parameter_2"):
                    for atlas_parameter_2 in parameter_2.iter("AtlasParameter2"):
                        if (
                            atlas_parameter_2.findall("Parameter_2_Value")[0].text
                            == inflation
                        ):
                            dwell_positions = (
                                atlas_parameter_2.findall("Dwell_Parameters")[0]
                                .findall("Dwell_Positions")[0]
                                .text.split(",")
                            )
                            dwell_times = (
                                atlas_parameter_2.findall("Dwell_Parameters")[0]
                                .findall("Dwell_Times")[0]
                                .text.split(",")
                            )
                            return [float(x) for x in dwell_positions], [
                                float(x) for x in dwell_times
                            ]


def geometry_function(distance: float):
    """Calculates geometry function for a single dwell position at a point of interest.

    Parameters
    ----------
    distance : float
        Distance from dwell position to point of interest.

    Returns
    -------
    float
        Geometry value accounting for distance from dwell position.

    Notes
    -----
    This currently assumes a point source, not a 2D source.
    """
    return 1 / math.pow(distance, 2)


def radial_function(distance: float):
    """Calculates radial dose function for a single dwell position at a point of interest.

    Parameters
    ----------
    distance : float
        Distance from dwell position to point of interest.

    Returns
    -------
    float
        Radial dose value accounting for attenuation of beam in water.
    """
    radial_data = pd.read_excel(
        os.path.join(_data_root_path, "xoft_axxent.xlsx"),
        sheet_name="Xoft Axxent",
        header=10,
        names=["r", "g"],
        usecols="B:C",
        nrows=14,
    )
    interpolant = interp1d(
        radial_data["r"].values, radial_data["g"].values, fill_value="extrapolate"
    )
    return interpolant(distance)


def anisotropy_function(distance: float, angle: float):
    """Calculates anisotropy function for a single dwell position at a point of interest.

    Parameters
    ----------
    distance : float
        Distance from dwell position to point of interest.

    Returns
    -------
    float
        Anisotropy value accounting for variations in dose with angle to point of interest.
    """
    anisotropy_data = pd.read_excel(
        os.path.join(_data_root_path, "xoft_axxent.xlsx"),
        sheet_name="Xoft Axxent",
        header=10,
        names=[
            "theta",
            "0.5",
            "1",
            "1.5",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "10",
            "12",
            "15",
        ],
        usecols="E:R",
        nrows=48,
    )
    return interpn(
        (
            anisotropy_data["theta"].values,
            [float(x) for x in anisotropy_data.columns[1:]],
        ),
        anisotropy_data.to_numpy()[:, 1:],
        (angle, distance),
    )[0]


def single_dwell_dose(distance: float, angle: float, time: float):
    """Calculate dose for a single dwell position at a point of interest.

    Parameters
    ----------
    distance : float
        Distance from dwell position to point of interest.
    angle : float
        Angle to point of interest from the distal catheter direction.
    time : float
        Dwell time (in sec).
        
    Returns
    -------
    float
        Dose at point (in Gy).

    """
    return (
        _reference_air_kerma
        * _dose_rate_constant
        * geometry_function(distance)
        * radial_function(distance)
        * anisotropy_function(distance, angle)
        / 60
        / 60
        / 100
        * time
    )


def plan_dose(
    plan_name: str,
    distance_from_balloon_surface: float = 0,
    angle_from_distal_direction: float = 90,
):
    """Calculate dose for a defined plan and point relative to balloon surface and distal direction.

    Parameters
    ----------
    plan_name : str
        Name of the template plan, for example '3-4_30cc'.
    distance_from_balloon_surface : float
        Distance to point of interest from the surface of the balloon (in cm, default 0).
    angle_from_distal_direction : float
        Angle to point of interest from the distal catheter direction (in degrees, default 90).

    Returns
    -------
    float
        Dose at point (in Gy).

    Notes
    -----
    The balloon radius is calculated assuming a spherical inflation, which is a source of uncertainty.
    """
    inflation = float(plan_name[4:-2])  # in cc
    dwell_positions, dwell_times = dwell_data(plan_name)
    # determine geometry relative to centre of assumed to be spherical balloon (approximation!)
    central_position = 25 - math.pow((3 * inflation) / (4 * math.pi), 1 / 3)
    distance_from_central_position = distance_from_balloon_surface + math.pow(
        (3 * inflation) / (4 * math.pi), 1 / 3
    )  # in cm
    offset_dwell_positions = (
        np.array(dwell_positions) - central_position
    )  # distal from centre is +ve, proximal is -ve, in cm
    dose_coordinates_from_central_position = (
        math.cos(math.radians(angle_from_distal_direction))
        * distance_from_central_position,
        math.sin(math.radians(angle_from_distal_direction))
        * distance_from_central_position,
    )
    # add doses for each dwell position
    cumulative_dose = 0
    for dwell in range(len(dwell_positions)):
        distance = math.sqrt(
            (dose_coordinates_from_central_position[0] - offset_dwell_positions[dwell])
            ** 2
            + dose_coordinates_from_central_position[1] ** 2
        )
        angle = math.degrees(
            math.acos(
                (
                    dose_coordinates_from_central_position[0]
                    - offset_dwell_positions[dwell]
                )
                / distance
            )
        )
        dose_from_dwell = (
            _reference_air_kerma
            * _dose_rate_constant
            * geometry_function(distance)
            * radial_function(distance)
            * anisotropy_function(distance, angle)
            / 60
            / 60
            / 100
            * dwell_times[dwell]
        )
        cumulative_dose += dose_from_dwell
    return cumulative_dose
