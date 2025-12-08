# -*- coding: utf-8 -*-
""" Kilovoltage dosimetry module.

This module provides kilovoltage dosimetry functions.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import modules
import os
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.interpolate import interpn

# find data location
_data_root_path = os.path.join(os.fspath(Path(__file__).parent.resolve()), "data")


def mass_absorption_coefficient(hvl: float, unit: str = "mm Al"):
    """Determine the ratios of average mass energy-absorption coefficients water to air.

    Parameters
    ----------
    hvl : float
        Half-value layer of the radiation.
    unit : str, optional
        Unit of the half-value layer, e.g. 'mm Al' (default) or 'mm Cu'.

    Returns
    -------
    float
        Ratio of average mass absorption coefficients water to air.
    """
    # pre-conditions
    if unit not in ["mm Al", "mm Cu"]:
        raise ValueError("Invalid filter specified.")
    if unit == "mm Al":
        hvl_values = np.load(os.path.join(_data_root_path, "tg61_table4a_hvl.npy"))
        mac_values = np.load(os.path.join(_data_root_path, "tg61_table4a_mac.npy"))
    else:
        hvl_values = np.load(os.path.join(_data_root_path, "tg61_table4b_hvl.npy"))
        mac_values = np.load(os.path.join(_data_root_path, "tg61_table4b_mac.npy"))
    interpolant = interp1d(hvl_values, mac_values, fill_value="extrapolate")
    return interpolant(hvl)


def backscatter_factor(
    hvl: float,
    unit: str = "mm Al",
    diameter: float = 5.0,
    ssd: float = 30.0,
    close_ended_correction: bool = False,
):
    """Determine the water-kerma based backscatter factor for a water phantom.

    Parameters
    ----------
    hvl : float
        Half-value layer of the radiation.
    unit : str, optional
        Unit of the half-value layer, e.g., 'mm Al' (default) or 'mm Cu'.
    diameter : float, optional
        Field diameter, in cm, e.g., 5 (default).
    ssd : float, optional
        Nominal focal source to surface distance, in cm, e.g., 30 (default).
    close_ended_correction : bool, optional
        Flag for multiplicative close-ended cone correction factor use (default is False).
    """
    # pre-conditions
    if unit not in ["mm Al", "mm Cu"]:
        raise ValueError("Invalid filter specified.")
    if unit == "mm Al" and close_ended_correction:
        raise ValueError("Close-ended cone corrections not possible for mm Al beams.")
    # For the data opened below, x is list of HVL values, y is list of diameters,
    # z is list of SSD values, and v is 3D array of data, as presented in TG-61
    if unit == "mm Al":
        if ssd < 10:
            x = np.load(os.path.join(_data_root_path, "tg61_table5a_x.npy"))
            y = np.load(os.path.join(_data_root_path, "tg61_table5a_y.npy"))
            z = np.load(os.path.join(_data_root_path, "tg61_table5a_z.npy"))
            v = np.load(os.path.join(_data_root_path, "tg61_table5a_v.npy"))
        else:
            x = np.load(os.path.join(_data_root_path, "tg61_table5b_x.npy"))
            y = np.load(os.path.join(_data_root_path, "tg61_table5b_y.npy"))
            z = np.load(os.path.join(_data_root_path, "tg61_table5b_z.npy"))
            v = np.load(os.path.join(_data_root_path, "tg61_table5b_v.npy"))
    else:
        x = np.load(os.path.join(_data_root_path, "tg61_table5c_x.npy"))
        y = np.load(os.path.join(_data_root_path, "tg61_table5c_y.npy"))
        z = np.load(os.path.join(_data_root_path, "tg61_table5c_z.npy"))
        v = np.load(os.path.join(_data_root_path, "tg61_table5c_v.npy"))
    factor = interpn((z, y, x), v, np.array([ssd, diameter, hvl]).T)
    if close_ended_correction:
        xi = np.load(os.path.join(_data_root_path, "tg61_table6_x.npy"))
        yi = np.load(os.path.join(_data_root_path, "tg61_table6_y.npy"))
        vi = np.load(os.path.join(_data_root_path, "tg61_table6_v.npy"))
        factor *= interpn((yi, xi), vi, np.array([diameter, hvl]).T)
    return factor
