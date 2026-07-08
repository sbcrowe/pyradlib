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
import os
import xml.etree.ElementTree as et
from functools import cached_property

import numpy as np
import polars as pl
import pydicom

# region RadixactPlan


class RadixactPlan:
    # region Constructors

    def __init__(self, ds: pydicom.Dataset) -> RadixactPlan:
        """Initialises the Radixact plan class.

        Parameters
        ----------
        ds : pydicom.Dataset
            DICOM RTPLAN dataset.

        Returns
        -------
        RadixactPlan
            The DICOM treatment plan wrapped in a helper class.
        """
        self._ds = ds

    @classmethod
    def from_dcm(cls, path: str | os.PathLike) -> RadixactPlan:
        """Reads treatment plan from DICOM RTPLAN file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the DICOM RTPLAN file.

        Returns
        -------
        RadixactPlan
            The DICOM treatment plan wrapped in a helper class.
        """
        ds = pydicom.dcmread(path)
        return cls(ds)

    # endregion

    # region Magic methods

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return f"RadixactPlan(SeriesInstanceUID={self._ds.SeriesInstanceUID})"

    # endregion

    # region Properties

    @cached_property
    def summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan parameters.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan parameters.
        """
        result = {}
        result["urn"] = self._ds.PatientID
        # TODO Manage names that aren't in format SMITH^JOHN
        result["last_name"] = str(self._ds.PatientName).split("^")[0]
        result["first_name"] = " ".join(str(self._ds.PatientName).split("^")[1:])
        result["physician"] = str(self._ds.OperatorsName)
        result["plan_name"] = self._ds.RTPlanName
        result["plan_date"] = self._ds.RTPlanDate
        result["plan_intent"] = self._ds.PlanIntent
        result["linear_accelerator"] = self._ds.BeamSequence[0].DeviceSerialNumber
        result["prescribed_dose"] = (
            float(self._ds.DoseReferenceSequence[0].TargetPrescriptionDose) * 100
        )
        result["number_of_fractions"] = int(
            self._ds.FractionGroupSequence[0].NumberOfFractionsPlanned
        )
        result["daily_dose"] = (
            float(self._ds.DoseReferenceSequence[0].TargetPrescriptionDose)
            * 100
            / int(self._ds.FractionGroupSequence[0].NumberOfFractionsPlanned)
        )
        result["beam_time"] = (
            float(
                self._ds.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset
            )
            * 60
        )
        result["pitch"] = float(self._ds.BeamSequence[0][0x300D, 0x1060].value)
        result["field_width"] = float(
            self._ds.BeamSequence[0]
            .BeamDescription.split(" ")[-1]
            .replace("width=", "")
        )
        result["dynamic_jaw"] = float(
            self._ds.BeamSequence[0]
            .ControlPointSequence[1]
            .BeamLimitingDevicePositionSequence[1]
            .LeafJawPositions[0]
        ) != float(
            self._ds.BeamSequence[0]
            .ControlPointSequence[-1]
            .BeamLimitingDevicePositionSequence[1]
            .LeafJawPositions[0]
        )
        result["gantry_period"] = float(self._ds.BeamSequence[0][0x300D, 0x1040].value)
        result["couch_speed"] = float(self._ds.BeamSequence[0][0x300D, 0x1080].value)
        result["control_points"] = int(self._ds.BeamSequence[0].NumberOfControlPoints)
        result["isocentre_iec_x"] = float(
            self._ds.BeamSequence[0].ControlPointSequence[1].IsocenterPosition[0]
        )
        result["isocentre_iec_y"] = float(
            self._ds.BeamSequence[0].ControlPointSequence[-1].IsocenterPosition[2]
        )
        result["isocentre_iec_z"] = float(
            self._ds.BeamSequence[0].ControlPointSequence[1].IsocenterPosition[1]
        )
        result["couch_travel"] = float(
            self._ds.BeamSequence[0].ControlPointSequence[1].IsocenterPosition[2]
        ) - float(
            self._ds.BeamSequence[0].ControlPointSequence[-1].IsocenterPosition[2]
        )
        return pl.DataFrame(result)

    # endregion


# endregion

# region RadixactPlanDetails


class RadixactPlanDetails:
    # region Constructors

    def __init__(self, df: pl.DataFrame) -> RadixactPlanDetails:
        """Initialise plan detail wrapper.

        Parameters
        ----------
        df : pl.DataFrame
            DataFrame containing plan details.

        Returns
        -------
        RadixactPlanDetails
            Plan details in helper wrapper.
        """
        self._df = df

    @classmethod
    def from_xml(cls, path: str | os.PathLike) -> RadixactPlanDetails:
        """Read plan details from an XML file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the plan detail xml file.

        Returns
        -------
        RadixactPlanDetails
            Plan details in helper wrapper.
        """
        details = {}
        tree = et.parse(path)
        root = tree.getroot()
        # Extract tomo plan details
        # fields = root.find("Model").find("Objects").find("Object").find("Fields")
        plan_details_objects = (
            root.find("Model")
            .find("Objects")
            .findall(".//*[@typeId='com.accuray.tps.domain_plan.TomoPlanDetails']")
        )
        if len(plan_details_objects) > 0:
            fields = plan_details_objects[0].find("Fields")
            details["machine_id"] = fields.find("TreatmentMachineId").text
            details["revision"] = fields.find("TreatmentMachineRevision").text
            details["type"] = fields.find("PlanDeliveryType").text
            details["delivery"] = fields.find("DeliveryScheme").text
            details["mode"] = fields.find("PlanMode").text
            details["couch_position"] = float(
                fields.find("CouchInsertionPositionMm").text
            )
            details["modulation_factor"] = float(
                fields.find("PlanningModulationFactor").text
            )
            details["pitch"] = float(fields.find("Pitch").text)
        # Extract imaging angles
        imaging_angle_objects = (
            root.find("Model")
            .find("Objects")
            .findall(
                ".//*[@typeId='com.accuray.tps.domain_plan.TomoMotionImagingAngle']"
            )
        )
        if len(imaging_angle_objects) > 0:
            imaging_angles = np.array(
                [
                    float(element.find("Fields").find("Angle").text)
                    for element in imaging_angle_objects
                ]
            )
            for index, angle in enumerate(imaging_angles):
                details[f"imaging_angle_{str(index + 1)}"] = angle
            for index in range(len(imaging_angles), 6, 1):
                details[f"imaging_angle_{str(index + 1)}"] = None
        # TODO Read other parameters, such as dose objectives and Synchrony imaging
        # angles.
        return cls(pl.DataFrame(details))

    # endregion


# endregion

# region RadixactPlanInformation


class RadixactPlanInformation:
    # region Constructors

    def __init__(
        self, df: pl.DataFrame, objectives: pl.DataFrame = None
    ) -> RadixactPlanInformation:
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
        self._df = df
        self._objectives = objectives

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
        profile_data["last_name"] = plan.find("PATIENT_PROFILE").find("LAST_NAME").text
        profile_data["first_name"] = (
            plan.find("PATIENT_PROFILE").find("FIRST_NAME").text
        )
        profile_data["plan_name"] = plan.find("PLAN_PROFILE").find("PLAN_NAME").text
        profile_data["date"] = (
            plan.find("PLAN_PROFILE").find("TIMESTAMP").find("DATETIME").text
        )
        profile_data["prescribed_dose"] = float(
            plan.find("PLAN_SETUP").find("PRESCRIBED_DOSE").text
        )
        profile_data["number_of_fractions"] = int(
            plan.find("PLAN_SETUP").find("NUMBER_OF_FRACTIONS").text
        )
        profile_data["number_of_vois"] = plan.find("VOI_INTERSECTION_SETTINGS").attrib[
            "size"
        ]
        profile_data["number_of_objectives"] = (
            plan.find("DX_VX_VALUES").find("DX_VX_DATA_SET").attrib["size"]
        )
        profile_data["number_of_fiducials"] = plan.find("FIDUCIALSET").attrib["size"]
        profile_data["ct_scanner"] = plan.find("DENSITY_MODEL").find("NAME").text
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
                metrics.append(objective_dict[spec_quantity] + spec_value)
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
        return cls(pl.DataFrame(profile_data), objectives_df)

    # endregion


# endregion


# region RadixactPlanSettings


class RadixactPlanSettings:
    # region Constructors

    def __init__(self, df: pl.DataFrame) -> RadixactPlanSettings:
        """Initialise plan settings class.

        Parameters
        ----------
        df : pl.DataFrame
            The plan settings dataframe.

        Returns
        -------
        RadixactPlanSettings
            Plan settings encapsulated in helping wrapper.
        """
        self._df = df

    @classmethod
    def from_plan_settings(cls, path: str | os.PathLike) -> RadixactPlanSettings:
        """Read plan settings from a plan.settings file.

        Parameters
        ----------
        path : path: str | os.PathLike
            Path to the plan.settings file.

        Returns
        -------
        RadixactPlanSettings
            Plan settings dataframe encapsulated in a RadixactPlanSettings wrapper.
        """
        data = {}
        with open(path, "r") as file:
            for line in file:
                if line.strip():
                    parameter, value = line.strip().split("=")
                    data[parameter.strip()] = value.strip()
        return cls(pl.DataFrame(data))

    # endregion

    # region Magic methods

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return f"RadixactPlanSettings(id={id(self)})"

    # endregion

    # region Properties

    @cached_property
    def summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan settings.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan settings.
        """
        result_dfs = []
        # Extract imaging protocols
        result_dfs.append(
            self._df.select(
                [
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(1)
                    .alias("kvct_protocol"),
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(2)
                    .alias("kvct_body_size"),
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(3)
                    .alias("kvct_field_of_view"),
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(4)
                    .alias("kvct_mode"),
                    (
                        pl.col("scanSelectionMaxEdge").cast(pl.Float64)
                        - pl.col("scanSelectionMinEdge").cast(pl.Float64)
                    ).alias("scan_length"),
                ]
            )
        )
        # Extract Synchrony data, if available angles
        if "measuredDifferenceTolerance" in self._df:
            result_dfs.append(
                self._df.select(
                    [
                        pl.col("measuredDifferenceTolerance")
                        .cast(pl.Float64)
                        .alias("measured_diff_tolerance"),
                        pl.col("potentialDifferenceTolerance")
                        .cast(pl.Float64)
                        .alias("potential_diff_tolerance"),
                        pl.col("targetOffsetTolerance")
                        .cast(pl.Float64)
                        .alias("target_offset_tolerance"),
                        pl.col("rigidBodyError")
                        .cast(pl.Float64)
                        .alias("rigid_body_tolerance"),
                        pl.col("trackingRange")
                        .cast(pl.Float64)
                        .alias("tracking_range"),
                        (
                            pl.col("TrackingViolationDelayMicroSeconds").cast(
                                pl.Float64
                            )
                            / 1000000
                        ).alias("auto_pause_delay"),
                        pl.col("sensitivity"),
                    ]
                )
            )
            result_dfs.append(
                self._df.select(
                    pl.col(imaging_angle_name)
                    .cast(pl.Float64)
                    .alias(f"imaging_angle_{str(index + 1)}")
                    if imaging_angle_name in self._df.columns
                    else pl.lit(None)
                    .cast(pl.Float64)
                    .alias(f"imaging_angle_{str(index + 1)}")
                    for index, imaging_angle_name in enumerate(
                        [
                            "RadiographAngle_0-Angle",
                            "RadiographAngle_1-Angle",
                            "RadiographAngle_2-Angle",
                            "RadiographAngle_3-Angle",
                            "RadiographAngle_4-Angle",
                            "RadiographAngle_5-Angle",
                        ]
                    )
                )
            )
        return pl.concat(result_dfs, how="horizontal")

    # endregion

    def to_csv(self, path: str | os.PathLike):
        """Wrie plan settings to comma-separated value file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the CSV file.
        """
        self._df.to_csv(path)


# endregion
