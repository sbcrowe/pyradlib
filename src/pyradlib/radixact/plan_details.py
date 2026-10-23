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
            The plan details as a DataFrame encapsulated with helper functions.
        """
        self._df: pl.DataFrame = df.fill_nan(None)

    @classmethod
    def from_csv(cls, path: str | os.PathLike) -> RadixactPlanDetails:
        """Read plan details from CSV file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the plan detail CSV file.

        Returns
        -------
        RadixactPlanDetails
            The plan details as a DataFrame encapsulated with helper functions.
        """
        return cls(pl.read_csv(path))

    @classmethod
    def from_path(cls, path: str | os.PathLike) -> RadixactPlanDetails:
        """Read plan details from path. Type of file identified by extension.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the plan detail file.

        Returns
        -------
        RadixactPlanDetails
            The plan details as a DataFrame encapsulated with helper functions.

        Raises
        ------
        ValueError
            If the extension is not "csv", "parquet", or "xml".
        """
        if path.endswith(".csv"):
            return cls.from_csv(path)
        elif path.endswith(".parquet"):
            return cls.from_parquet(path)
        elif path.endswith(".xml"):
            return cls.from_xml(path)
        else:
            raise ValueError(
                f"Unsupported extension for RadixactPlanDetails: {path.split('.')[-1]}"
            )

    @classmethod
    def from_parquet(cls, path: str | os.PathLike) -> RadixactPlanDetails:
        """Read plan details from Parquet file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the plan detail Parquet file.

        Returns
        -------
        RadixactPlanDetails
            The plan details as a DataFrame encapsulated with helper functions.
        """
        return cls(pl.read_parquet(path))

    @classmethod
    def from_xml(cls, path: str | os.PathLike) -> RadixactPlanDetails:
        """Read plan details from XML file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the plan detail XML file.

        Returns
        -------
        RadixactPlanDetails
            The plan details as a DataFrame encapsulated with helper functions.
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

    # region Magic methods

    def __len__(self) -> int:
        """Returns number of plan details.

        Returns
        -------
        int
            Number of plan details.
        """
        return len(self._df)

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return self._df.__repr__()

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

    # region Public methods

    def to_csv(self, path: str | os.PathLike) -> None:
        """Writes plan details to CSV file.

        Parameters
        ----------
        path : str | os.PathLike
            Path for CSV file to be written.
        """
        self._df.write_csv(path)

    def to_parquet(self, path: str | os.PathLike) -> None:
        """Writes plan details to Parquet file.

        Parameters
        ----------
        path : str | os.PathLike
            Path for Parquet file to be written.
        """
        self._df.write_parquet(path)

    def to_path(self, path: str | os.PathLike) -> None:
        """Writes plan details to file. Type of file identified by extension.

        Parameters
        ----------
        path : str | os.PathLike
            Path for file to be written.
        """
        if path.endswith(".csv"):
            self.to_csv(path)
        elif path.endswith(".parquet"):
            self.to_parquet(path)
        else:
            raise ValueError(
                f"Unsupported extension for RadixactTiming: {path.split('.')[-1]}"
            )

    # endregion


# endregion
