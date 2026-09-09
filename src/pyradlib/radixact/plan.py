"""Plan analysis module.

This module provides functionality for processing of general planning data from
Radixact treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import os
from functools import cached_property

import polars as pl
import pydicom


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
        self._ds: pydicom.Dataset = ds

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
        if ds.Modality != "RTPLAN":
            raise ValueError("'Dataset' object attribute 'Modality' is not 'RTPLAN'")
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
        result["patient_name"] = str(self._ds.PatientName)
        result["physician_name"] = str(self._ds.OperatorsName)
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

    # region Public methods

    # endregion
