import os

from pyradlib.radixact.dataset import RadixactDataset
from pyradlib.radixact.motion import RadixactSynchronyMotion
from pyradlib.radixact.plan import (
    RadixactPlan,
    RadixactPlanDetails,
    RadixactPlanInformation,
    RadixactPlanSettings,
)
from pyradlib.radixact.sinogram import RadixactSinogram
from pyradlib.radixact.timing import RadixactTiming

__all__ = [
    "RadixactDataset",
    "RadixactPlan",
    "RadixactPlanDetails",
    "RadixactPlanInformation",
    "RadixactPlanSettings",
    "RadixactSinogram",
    "RadixactSynchronyMotion",
    "RadixactTiming",
]


def read_dir(path: str):
    if os.path.exists(os.path.join(path, "PatientExportDataBO.xml")):
        return RadixactDataset.from_data_extractor(path)
    else:
        return RadixactDataset.from_delivery_analysis(path)


def read_da(path: str):
    return RadixactDataset.from_delivery_analysis(path)


def read_dat(path: str):
    return RadixactTiming.from_dat(path)


def read_det(path: str):
    """Reads detector sinogram from a .det file.

    Parameters
    ----------
    path : str
        Path to the .det file.
    clip : bool, optional
        Flag to indicate whether non-detector channels should be removed.
    flip : bool, optional
        Flag to indicate whether channels should be mirrored left-to-right,
        for ease of visual comparison against leaf sinograms.

    Returns
    -------
    Sinogram
        The det file is returned as a two-dimensional numpy array encapsulated
        with helper functions.

    Notes
    -----
    The function converts the binary-to-ASCII (hexadecimal) encoded data to a
    sinogram. This depends on the presence of fragment counts to determine how
    many projections exist for each leaf channel. For each leaf channel, the
    number of projections with open leaves is calculated from paired '>f8' tau
    values.
    """
    return RadixactSinogram.from_det(path)


def read_dplan(path: str):
    """Reads planned or telemetry sinogram from a .dplan file.

    Parameters
    ----------
    path : str
        Path to the .dplan file.

    Returns
    -------
    Sinogram
        The det file is returned as a two-dimensional numpy array encapsulated
        with helper functions.

    Notes
    -----
    The function converts the binary-to-ASCII (hexadecimal) encoded data to a
    sinogram. This depends on the presence of fragment counts to determine how
    many projections exist for each leaf channel. For each leaf channel, the
    number of projections with open leaves is calculated from paired '>f8' tau
    values.
    """
    return RadixactSinogram.from_dplan(path)


def read_motion(path: str):
    return RadixactSynchronyMotion.from_xml(path)


def read_pde(path: str):
    return RadixactDataset.from_data_extractor(path)


def read_plan(path: str):
    return RadixactPlan.from_dcm(path)


def read_plan_details(path: str):
    return RadixactPlanDetails.from_xml(path)


def read_plan_information(path: str):
    return RadixactPlanInformation.from_xml(path)


def read_plan_settings(path: str):
    return RadixactPlanSettings.from_plan_settings(path)
