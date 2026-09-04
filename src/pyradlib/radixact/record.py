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

import matplotlib as mpl
import polars as pl
import pydicom
import seaborn as sns


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
        self._ds: pydicom.Dataset = ds

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
        if ds.Modality != "RTRECORD":
            raise ValueError("'Dataset' object attribute 'Modality' is not 'RTRECORD'")
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
        result["linear_accelerator"] = self._ds.TreatmentMachineSequence[
            0
        ].DeviceSerialNumber
        # if session has beam sequence
        if "TreatmentSessionBeamSequence" in self._ds:
            result["fraction"] = self._ds.TreatmentSessionBeamSequence[
                0
            ].CurrentFractionNumber
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
            # if session has corrections
            description = self._ds.TreatmentSessionBeamSequence[
                0
            ].BeamDescription.split("\n")
            if len(description) > 6:
                result["correction_x"] = float(description[1].split("=")[1])
                result["correction_y"] = float(description[2].split("=")[1])
                result["correction_z"] = float(description[3].split("=")[1])
                result["correction_roll"] = float(description[5].split("=")[1])
            else:
                result["correction_x"] = None
                result["correction_y"] = None
                result["correction_z"] = None
                result["correction_roll"] = None
        else:
            result["fraction"] = None
            result["planned_beam_time"] = None
            result["delivered_beam_time"] = None
            result["planned_mu"] = None
            result["delivered_mu"] = None
            result["correction_x"] = None
            result["correction_y"] = None
            result["correction_z"] = None
            result["correction_roll"] = None
        return pl.DataFrame(result)

    # endregion

    # region Static methods

    @staticmethod
    def plot_correction_histogram(
        summary: pl.DataFrame,
        parameters: list[str] | None = None,
        binwidth: float = 0.25,
        col_wrap: int = 4,
        sharex: bool = False,
        sharey: bool = False,
    ) -> mpl.Figure:
        """Generate a histogram plot of treatment couch adjustments after registration.

        Parameters
        ----------
        summary : pl.DataFrame
            The summary of delivery records.
        parameters : list[str], optional
            The couch adjustments to include in the figure. Possible options include:
            "correction_x", "correction_y", "correction_z"
                Couch adjustment in the IEC-X, IEC-Y and/or IEC-Z dimensions.
            "correction_roll"
                Gantry angle adjustment.
            Default is None, in which case correction_x, correction_y, correction_z,
            and correction_roll are used.
        binwidth: float, optional.
            Width of bin, in mm or degrees. Default is 0.25 (mm or degrees).
        col_wrap : int, optional
            Number of column facets within a row. Default is 4.
        sharex : bool, optional
            Flag for whether facets should share x axes. Default is False.
        sharey : bool, optional
            Flag for whether facets should share y axes. Default is False.

        Returns
        -------
        mpl.Figure
            Histogram of target offset values in each dimension.
        """
        if parameters is None:
            parameters = [
                "correction_x",
                "correction_y",
                "correction_z",
                "correction_roll",
            ]
        mapping = {
            "correction_x": "IEC-X",
            "correction_y": "IEC-Y",
            "correction_z": "IEC-Z",
            "correction_roll": "Roll",
        }
        indices = [
            index for index in ["session_index", "patient_index"] if index in summary
        ]
        select = [pl.col(col) for col in indices + parameters if col in summary]
        unpivot_df = (
            summary.select(select)
            .unpivot(index=indices)
            .select(
                [
                    (pl.col("session_index") + 1).alias("Session"),
                    pl.col("variable").replace(mapping).alias("Parameter"),
                    pl.col("value").alias("Correction (mm or °)"),
                ]
            )
        )
        ax = sns.FacetGrid(
            unpivot_df, col="Parameter", col_wrap=col_wrap, sharex=sharex, sharey=sharey
        )
        ax.map(sns.histplot, "Correction (mm or °)", binwidth=binwidth)
        ax.set_titles(template="{col_name}")
        for a in ax.axes.flat:
            a.set_xlabel("Correction (mm)")
        ax.axes.flat[-1].set_xlabel("Correction (°)")
        return ax.figure

    # endregion


# endregion
