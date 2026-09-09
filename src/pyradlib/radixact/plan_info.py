"""Plan information module.

This module provides functionality for processing of plan information data from
Radixact treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import datetime
import os
import xml.etree.ElementTree as et

import polars as pl


class RadixactPlanInformation:
    # region Constructors

    def __init__(self, df: pl.DataFrame) -> RadixactPlanInformation:
        """Initialises a plan information wrapper.

        Parameters
        ----------
        df : pl.DataFrame
            Plan information dataframe.
        objectives : pl.DataFrame, optional
            Dose objectives dataframe. Default is default None.

        Returns
        -------
        RadixactPlanInformation
            Plan information data, in helper wrapper.
        """
        self._df: pl.DataFrame = df

    @classmethod
    def from_xml(cls, path: str | os.PathLike) -> RadixactPlanInformation:
        """Reads general plan and dose objective information from a GeneralPlan*.xml
        file.

        Parameters
        ----------
        path: str | os.PathLike
            Path to the GeneralPlan*.xml file.

        Returns
        -------
        RadixactPlanInformation
            Plan profile and dose objective information encapsulated in a
            RadixactPlanInformation wrapper.
        """
        profile_data = {}
        tree = et.parse(path)
        root = tree.getroot()
        plan = root.find("GENERAL_PLAN")
        profile_data["urn"] = plan.find("PATIENT_PROFILE").find("MEDICAL_ID").text
        profile_data["patient_name"] = (
            plan.find("PATIENT_PROFILE").find("LAST_NAME").text
            + "^"
            + plan.find("PATIENT_PROFILE").find("FIRST_NAME").text
        )
        profile_data["plan_name"] = plan.find("PLAN_PROFILE").find("PLAN_NAME").text
        profile_data["datetime"] = datetime.datetime.fromisoformat(
            plan.find("PLAN_PROFILE")
            .find("TIMESTAMP")
            .find("DATETIME")
            .text.replace("TZ", "")
        )
        profile_data["prescribed_dose"] = float(
            plan.find("PLAN_SETUP").find("PRESCRIBED_DOSE").text
        )
        profile_data["number_of_fractions"] = int(
            plan.find("PLAN_SETUP").find("NUMBER_OF_FRACTIONS").text
        )
        profile_data["number_of_vois"] = int(
            plan.find("VOI_INTERSECTION_SETTINGS").attrib["size"]
        )
        profile_data["number_of_objectives"] = int(
            plan.find("DX_VX_VALUES").find("DX_VX_DATA_SET").attrib["size"]
        )
        profile_data["number_of_fiducials"] = int(
            plan.find("FIDUCIALSET").attrib["size"]
        )
        profile_data["ct_scanner"] = plan.find("DENSITY_MODEL").find("NAME").text
        return cls(pl.DataFrame(profile_data))

    # endregion

    # region Properties

    @property
    def summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan information.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan information.
        """
        return self._df

    # endregion


# region RadixactDoseObjectives


class RadixactDoseObjectives:
    # region Constructors

    def __init__(self, df: pl.DataFrame) -> RadixactDoseObjectives:
        """Initialises a planning dose objectives wrapper.

        Parameters
        ----------
        df : pl.DataFrame
            Dose objectives dataframe.

        Returns
        -------
        RadixactDoseObjectives
            Dose objective data, in helper wrapper.
        """
        self._df: pl.DataFrame = df

    @classmethod
    def from_xml(cls, path: str | os.PathLike) -> RadixactDoseObjectives:
        tree = et.parse(path)
        root = tree.getroot()
        plan = root.find("GENERAL_PLAN")

        # Read dose objectives
        vois = []
        metrics = []
        objectives = []
        voi_dict = {}
        for voi in root.find("GENERAL_PLAN").find("AUTOSEG").find("VOISET"):
            voi_dict[voi.attrib["id"]] = voi.attrib["name"]
        unit_dict = {"1": "cGy", "2": "%", "4": "%"}
        objective_dict = {"1": "V", "2": "V", "4": "D"}
        operator_dict = {"1": "<", "3": ">"}
        for dx_vx in (
            root.find("GENERAL_PLAN").find("DX_VX_VALUES").find("DX_VX_DATA_SET")
        ):
            voi = voi_dict[dx_vx.find("ACTIVE_PLAN_VOI_INDEX").text]
            spec_quantity = dx_vx.find("SPECIFIED_QUANTITY").text
            spec_value = round(float(dx_vx.find("SPECIFIED_VALUE").text))
            if dx_vx.find("DX_VX_CRITERIA") is None:
                vois.append(voi)
                metrics.append(objective_dict[spec_quantity] + str(spec_value))
                objectives.append("undefined")
            else:
                crit_quantity = (
                    dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_QUANTITY").text
                )
                crit_operator = dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_OP").text
                crit_value1 = (
                    dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_OPERAND1").text
                )
                crit_value2 = (
                    dx_vx.find("DX_VX_CRITERIA").find("CRITERIA_OPERAND2").text
                )
                if crit_operator == "5":
                    vois.append(voi)
                    metrics.append(
                        objective_dict[spec_quantity]
                        + str(spec_value)
                        + unit_dict[spec_quantity]
                    )
                    objectives.append("> " + crit_value1 + unit_dict[crit_quantity])
                    vois.append(voi)
                    metrics.append(
                        objective_dict[spec_quantity]
                        + str(spec_value)
                        + unit_dict[spec_quantity]
                    )
                    objectives.append("< " + crit_value2 + unit_dict[crit_quantity])
                else:
                    vois.append(voi)
                    metrics.append(
                        objective_dict[spec_quantity]
                        + str(spec_value)
                        + unit_dict[spec_quantity]
                    )
                    objectives.append(
                        operator_dict[crit_operator]
                        + " "
                        + crit_value1
                        + unit_dict[crit_quantity]
                    )
        objectives_df = pl.DataFrame(
            {"voi": vois, "metric": metrics, "objective": objectives}
        )
        return cls(objectives_df)

    # endregion


# endregion
