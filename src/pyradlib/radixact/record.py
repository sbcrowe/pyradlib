"""Record module.

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
from functools import cached_property

import polars as pl
import pydicom


class RadixactRecord:
    # region Constructors

    def __init__(self, ds: pydicom.Dataset) -> RadixactRecord:
        """Initialises the Radixact plan class.

        Parameters
        ----------
        ds : pydicom.Dataset
            DICOM RTRECORD dataset.

        Returns
        -------
        RadixactRecord
            The DICOM treatment plan wrapped in a helper class.
        """
        self._ds = ds

    @classmethod
    def from_dcm(cls, path: str | os.PathLike) -> RadixactRecord:
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
        return f"RadixactRecord(SeriesInstanceUID={self._ds.SeriesInstanceUID})"

    # endregion

    # region Properties

    @cached_property
    def summary(self) -> pl.DataFrame:
        """Produce summary of treatment record parameters.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment record parameters.
        """
        result = {}
        result["urn"] = self._ds.PatientID
        # TODO Manage names that aren't in format SMITH^JOHN
        result["last_name"] = str(self._ds.PatientName).split("^")[0]
        result["first_name"] = " ".join(str(self._ds.PatientName).split("^")[1:])
        result["operator"] = str(self._ds.OperatorsName)
        # result["plan_name"] = self._ds.RTPlanName
        result["treatment_date"] = self._ds.TreatmentDate
        result["treatment_time"] = self._ds.TreatmentTime
        # result["plan_intent"] = self._ds.PlanIntent
        result["plan_uid"] = self._ds.ReferencedRTPlanSequence[
            0
        ].ReferencedSOPInstanceUID
        result["fraction"] = self._ds.TreatmentSessionBeamSequence[
            0
        ].CurrentFractionNumber
        result["linear_accelerator"] = self._ds.TreatmentMachineSequence[
            0
        ].DeviceSerialNumber
        result["planned_beam_time"] = self._ds.TreatmentSessionBeamSequence[
            0
        ].SpecifiedTreatmentTime
        result["delivered_beam_time"] = self._ds.TreatmentSessionBeamSequence[
            0
        ].DeliveredTreatmentTime
        result["planned_mu"] = self._ds.TreatmentSessionBeamSequence[
            0
        ].SpecifiedSecondaryMeterset
        result["delivered_mu"] = self._ds.TreatmentSessionBeamSequence[
            0
        ].DeliveredSecondaryMeterset
        description = self._ds.TreatmentSessionBeamSequence[0].BeamDescription.split(
            "\n"
        )
        result["couch_adjustment_x"] = float(description[1].split("=")[1])
        result["couch_adjustment_y"] = float(description[2].split("=")[1])
        result["couch_adjustment_z"] = float(description[3].split("=")[1])
        result["couch_adjustment_roll"] = float(description[5].split("=")[1])
        return pl.DataFrame(result)

    # endregion


# endregion
