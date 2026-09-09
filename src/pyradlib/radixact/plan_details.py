"""Plan details module.

This module provides functionality for processing of plan details data from
Radixact treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import os
import xml.etree.ElementTree as et

import numpy as np
import polars as pl


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
        self._df : pl.DataFrame = df

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
            details["revision"] = int(fields.find("TreatmentMachineRevision").text)
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
                details[f"imaging_angle_{index + 1}"] = angle
            for index in range(len(imaging_angles), 6, 1):
                details[f"imaging_angle_{index + 1}"] = None
        # TODO Read other parameters, such as dose objectives and Synchrony imaging
        # angles.
        return cls(pl.DataFrame(details))

    # endregion

    # region Properties

    @property
    def summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan details.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan details.
        """
        return self._df

    # endregion


# endregion
