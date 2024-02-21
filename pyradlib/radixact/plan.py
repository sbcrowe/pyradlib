# -*- coding: utf-8 -*-
"""Plan analysis module.

This module provides functionality for processing of general planning data from 
Synchrony treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import numpy as np
import numpy.typing as npt
import pandas as pd
import xml.etree.ElementTree as et


def read_plan_data(xml_filepaths: npt.ArrayLike):
    """Extracts data from multiple GeneralPlan*.xml files produced for the Radixact
    Delivery Analysis tool.

    Parameters
    ----------
    xml_filepaths : array_like
        List of paths to XML files containing general plan data.

    Returns
    -------
    results : Pandas DataFrame
        A dataframe containing data extracted from the xml files.

    Notes
    -----
    These files are cached in C:/tomo/da/pts/URnumber/*motionData.xml when patient data
    is loaded within the Delivery Analysis tool. There may be multiple files per fraction.
    """
    header = [
        "URN",
        "Last Name",
        "First Name",
        "Plan",
        "Date",
        "Prescribed dose (cGy)",
        "Number of fractions",
        "Number of VOIs",
        "Number of optimisation objectives",
        "Number of fiducials",
        "CT scanner",
    ]
    results = []
    for xml_path in xml_filepaths:
        curr_result = []
        tree = et.parse(xml_path)
        root = tree.getroot()
        plan = root.find("GENERAL_PLAN")
        curr_result.append(plan.find("PATIENT_PROFILE").find("MEDICAL_ID").text)
        curr_result.append(plan.find("PATIENT_PROFILE").find("LAST_NAME").text)
        curr_result.append(plan.find("PATIENT_PROFILE").find("FIRST_NAME").text)
        curr_result.append(plan.find("PLAN_PROFILE").find("PLAN_NAME").text)
        curr_result.append(
            plan.find("PLAN_PROFILE").find("TIMESTAMP").find("DATETIME").text
        )
        curr_result.append(float(plan.find("PLAN_SETUP").find("PRESCRIBED_DOSE").text))
        curr_result.append(
            int(plan.find("PLAN_SETUP").find("NUMBER_OF_FRACTIONS").text)
        )
        curr_result.append(plan.find("VOI_INTERSECTION_SETTINGS").attrib["size"])
        curr_result.append(
            plan.find("DX_VX_VALUES").find("DX_VX_DATA_SET").attrib["size"]
        )
        curr_result.append(plan.find("FIDUCIALSET").attrib["size"])
        curr_result.append(plan.find("DENSITY_MODEL").find("NAME").text)
        results.append(curr_result)
    return pd.DataFrame(results, columns=header)


def read_dose_objectives(xml_filepath: str):
    """Extracts dose objective data from a GeneralPlan*.xml file produced for the
    Radixact Delivery Analysis tool.

    Parameters
    ----------
    xml_filepath : str
        Path to XML files containing general plan data.

    Returns
    -------
    results : array_like
        A list containing dose objectives.

    Notes
    -----
    These files are cached in C:/tomo/da/pts/URnumber/*motionData.xml when patient data
    is loaded within the Delivery Analysis tool. There may be multiple files per fraction.
    """
    results = []
    tree = et.parse(xml_filepath)
    root = tree.getroot()
    voi_dict = {}
    for voi in root.find("GENERAL_PLAN").find("AUTOSEG").find("VOISET"):
        voi_dict[voi.attrib["id"]] = voi.attrib["name"]
    unit_dict = {"1": "cGy", "2": "%", "4": "%"}
    objective_dict = {"1": "V", "2": "V", "4": "D"}
    operator_dict = {"1": "<", "3": ">"}
    for dx_vx in root.find("GENERAL_PLAN").find("DX_VX_VALUES").find("DX_VX_DATA_SET"):
        voi = voi_dict[dx_vx.find("ACTIVE_PLAN_VOI_INDEX").text]
        spec_quantity = dx_vx.find("SPECIFIED_QUANTITY").text
        spec_value = dx_vx.find("SPECIFIED_VALUE").text
        if dx_vx.find("DX_VX_CRITERIA") is None:
            results.append(
                voi
                + " "
                + objective_dict[spec_quantity]
                + spec_value
                + unit_dict[spec_quantity]
                + " is undefined"
            )
        else:
            crit_quantity = dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_QUANTITY").text
            crit_operator = dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_OP").text
            crit_value1 = dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_OPERAND1").text
            crit_value2 = dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_OPERAND2").text
            if crit_operator == 5:
                results.append(
                    voi
                    + " "
                    + objective_dict[spec_quantity]
                    + spec_value
                    + unit_dict[spec_quantity]
                    + " in ["
                    + crit_value1
                    + unit_dict[crit_quantity]
                    + ", "
                    + crit_value2
                    + unit_dict[crit_quantity]
                    + "]"
                )
            else:
                results.append(
                    voi
                    + " "
                    + objective_dict[spec_quantity]
                    + spec_value
                    + unit_dict[spec_quantity]
                    + " "
                    + operator_dict[crit_operator]
                    + " "
                    + crit_value1
                    + unit_dict[crit_quantity]
                )
    return results
