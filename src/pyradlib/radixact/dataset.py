"""RadixactDataset module.

This module provides functionality for processing of data from Radixact treatments.
"""

from __future__ import annotations

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import glob
import logging
import os
import shutil
import xml.etree.ElementTree as et
from functools import cached_property
from typing import ClassVar

import matplotlib as mpl
import numpy as np
import polars as pl
import pydicom
import seaborn as sns

from pyradlib.radixact.motion import RadixactSynchronyMotion
from pyradlib.radixact.plan import (
    RadixactPlan,
    RadixactPlanDetails,
    RadixactPlanInformation,
    RadixactPlanSettings,
)
from pyradlib.radixact.sinogram import RadixactSinogram
from pyradlib.radixact.timing import RadixactTiming

logger = logging.getLogger(__name__)


class RadixactDataset:
    # region Class variables

    _COMPRESSED_EXTENSIONS: ClassVar[dict[str, str]] = {
        "Detector Sinogram": "parquet",
        "Motion Data": "npz",
        "Plan Sinogram": "parquet",
        "Telemetry Sinogram": "parquet",
    }

    _DEFAULT_EXTENSIONS: ClassVar[dict[str, str]] = {
        "Detector Sinogram": "det",
        "Dose Distribution": "dcm",
        "Motion Data": "xml",
        "Plan": "dcm",
        "Plan Details": "xml",
        "Plan Information": "xml",
        "Plan Settings": "txt",
        "Plan Sinogram": "dplan",
        "Record": "dcm",
        "Structure Set": "dcm",
        "Telemetry Sinogram": "dplan",
        "Telemetry Timing": "dat",
    }

    _DEFAULT_PREFIXES: ClassVar[dict[str, str]] = {
        "Detector Sinogram": "detectorfluence",
        "Dose Distribution": "rtdose",
        "Motion Data": "motiondata",
        "Plan": "rtplan",
        "Plan Details": "plandetails",
        "Plan Information": "planinfo",
        "Plan Settings": "plansettings",
        "Plan Sinogram": "planfluence",
        "Record": "rtrecord",
        "Structure Set": "rtss",
        "Telemetry Sinogram": "telemfluence",
        "Telemetry Timing": "telemtiming",
    }

    # endregion

    # region Constructors

    def __init__(self, files: pl.DataFrame) -> RadixactDataset:
        """Initialises an object corresponding to a Radixact patient dataset.

        Parameters
        ----------
        files : pl.DataFrame
            DataFrame containing the files comprising the dataset.

        Returns
        -------
        RadixactDataset
            The encapsulated dataset object.
        """
        self._files = files
        logger.info(f"Dataset initialised with {len(self._files)!s} files")

    @classmethod
    def from_data_extractor(cls, path: str | os.PathLike) -> RadixactDataset:
        """Reads dataset from a Patient Data Extractor folder.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the patient data folder.

        Returns
        -------
        RadixactDataset
            The encapsulated dataset object.

        Notes
        -----
        The specified path should correspond to a Patient Data Extractor folder,
        containing a PatientExportDataBO.xml file.
        """
        dataset = cls(cls._extractor_data_files_dataframe(path))
        return dataset

    @classmethod
    def from_delivery_analysis(cls, path: str | os.PathLike) -> RadixactDataset:
        """Reads dataset from a Delivery Analysis folder.

        Parameters
        ----------
        path : str | os.PathLike
            Path containing the Delivery Analysis cache files to be read.

        Returns
        -------
        RadixactDataset
            The encapsulated dataset object.

        Notes
        -----
        The specified directory should correspond to a Delivery Analysis cache
        folder, e.g., containing DetectorByProj*.det, TelemFluence_*.dplan, and
        FinalDeliveryPlan*.dplan files.
        """
        return cls.from_path_list(
            sorted([os.path.join(path, filename) for filename in os.listdir(path)])
        )

    @classmethod
    def from_path(cls, path: str | os.PathLike) -> RadixactDataset:
        """Reads dataset from a directory.

        Parameters
        ----------
        path : str | os.PathLike
            Path of directory containing the files to be read.

        Returns
        -------
        RadixactDataset
            The encapsulated dataset object.

        Notes
        -----
        The method will check whether the directory is an anonymised patient data
        extractor export, and if not, assume it is a delivery analysis directory.
        """
        if os.path.exists(os.path.join(path, "PatientExportDataBO.xml")):
            return RadixactDataset.from_data_extractor(path)
        else:
            return RadixactDataset.from_delivery_analysis(path)

    @classmethod
    def from_path_list(cls, paths: list[str] | list[os.PathLike]) -> RadixactDataset:
        """
        Parameters
        ----------
        paths : list[str] | list[os.PathLike]
            Paths of the files to be included in the DataSet.

        Returns
        -------
        RadixactDataset
            The encapsulated dataset object.
        """
        return cls(
            cls._categorise_data_files_dataframe(
                pl.DataFrame({"curr_path": paths}), series_name="curr_path"
            )
        )

    # endregion

    # region Magic methods

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return f"RadixactDataset(num_files={len(self._files)})"

    # endregion

    # region Properties

    @property
    def combined_plan_summary(self) -> pl.DataFrame:
        """Return joined plan, plan detail and plan information summaries.

        Returns
        -------
        pl.DataFrame
            Combined plan, plan detail and plan information summaries.
        """
        shared_plan_detail_columns = list(
            set(self.plan_summary.columns) & set(self.plan_details_summary.columns)
        )
        shared_plan_info_columns = list(
            set(self.plan_summary.columns) & set(self.plan_informations_summary.columns)
        )
        return self.plan_summary.join(
            self.plan_details_summary, on=shared_plan_detail_columns, how="inner"
        ).join(self.plan_informations[0]._df, on=shared_plan_info_columns, how="inner")

    @cached_property
    def detector_sinograms(self) -> list[RadixactSinogram]:
        """Returns list containing detector sinograms.

        Returns
        -------
        list[RadixactSinogram]
            List of detector sinograms for each delivery session in the dataset.
        """
        # TODO Potentially include orig_path or curr_path.split("-")[-3] as plan_uid,
        # [-2] as fraction number, and [-1].replace(".det", "") as fragment.
        detector_sinograms = []
        detector_sinogram_files = self._get_filtered_series_list(
            "Detector Sinogram", series="curr_path"
        )
        logger.debug(
            f"Dataset contains {len(detector_sinogram_files)} detector sinogram files"
        )
        if len(detector_sinogram_files) == 0:
            logger.debug("Dataset contains no detector sinogram files")
            return detector_sinograms
        else:
            for detector_sinogram_file in detector_sinogram_files:
                logger.debug(f"Loading {detector_sinogram_file}")
                detector_sinograms.append(
                    RadixactSinogram.from_det(detector_sinogram_file)
                )
            logger.info(f"Loaded {len(detector_sinograms)} detector sinogram files")
            return detector_sinograms

    @cached_property
    def dose_distributions(self) -> list[pydicom.Dataset]:
        """Returns list containing dose distributions.

        Returns
        -------
        list[pydicom.Dataset]
            List of dose distributions for each plan in the dataset.
        """
        # TODO Potentially include orig_path or curr_path.split("-")[-2] as fraction number
        # or include entire UID from orig_path or curr_path.
        dose_distributions = []
        dose_distribution_files = self._get_filtered_series_list(
            "Dose Distribution", series="curr_path"
        )
        logger.debug(
            f"Dataset contains {len(dose_distribution_files)} dose distribution files"
        )
        if len(dose_distribution_files) == 0:
            logger.debug("Dataset contains no dose distribution files")
            return dose_distributions
        else:
            for dose_distribution_file in dose_distribution_files:
                logger.debug(f"Loading {dose_distribution_file}")
                dose_distributions.append(pydicom.dcmread(dose_distribution_file))
            logger.info(f"Loaded {len(dose_distributions)} dose distribution files")
            return dose_distributions

    @cached_property
    def motion(self) -> RadixactSynchronyMotion:
        """Returns concatenated motion data.

        Returns
        -------
        RadixactSynchronyMotion
            Concated Synchrony motion data from the dataset.
        """
        if len(self.motions) == 0:
            return None
        else:
            return RadixactSynchronyMotion.from_session_motions(self.motions)

    @cached_property
    def motion_metrics(self) -> pl.DataFrame:
        return self.motion.metrics

    @cached_property
    def motions(self) -> list[RadixactSynchronyMotion]:
        """Returns list containing motion data.

        Returns
        -------
        list[RadixactSynchronyMotion]
            List of Synchrony motion data from each delivery session in the dataset.
        """
        # TODO Potentially include orig_path or curr_path.split("-")[-2] as fraction number
        # or include entire UID from orig_path or curr_path.
        motions = []
        motion_files = self._get_filtered_series_list("Motion Data", series="curr_path")
        if len(motion_files) == 0:
            logger.debug("Dataset contains no motion data files")
            return motions
        else:
            for motion_file in motion_files:
                logger.debug(f"Loading {motion_file}")
                motion = RadixactSynchronyMotion.from_xml(motion_file)
                if motion is not None:
                    motions.append(motion)
            logger.info(f"Loaded {len(motions)} motion data files")
            return motions

    @cached_property
    def plan_details_summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan details, for all plans in the dataset.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan details, for all plans in the
            dataset.
        """
        return pl.concat(
            [
                plan_details.summary.with_columns(plan_index=pl.lit(key))
                for key, plan_details in enumerate(self.plan_details)
            ]
        ).select([pl.col("plan_index"), pl.all().exclude("plan_index")])

    @cached_property
    def plan_informations_summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan information, for all plans in the dataset.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan information, for all plans in the
            dataset.
        """
        return pl.concat(
            [
                plan_information.summary.with_columns(plan_set_index=pl.lit(key))
                for key, plan_information in enumerate(self.plan_informations)
            ]
        ).select([pl.col("plan_set_index"), pl.all().exclude("plan_set_index")])

    @cached_property
    def plans(self) -> list[RadixactPlan]:
        """Reurns list containing treatment plans.

        Returns
        -------
        list[RadixactPlan]
            List of treatment plans in the dataset.
        """
        # TODO Potentially extract ds.BeamSequence[0].DeviceSerialNumber and
        # ds.RTPlanLabel to be used as identifiers in a dictionary.
        plans = []
        plan_files = self._get_filtered_series_list("Plan", series="curr_path")
        if len(plan_files) == 0:
            logger.debug("Dataset contains no plan files")
            return plans
        else:
            for plan_file in plan_files:
                logger.debug(f"Loading {plan_file}")
                plans.append(RadixactPlan.from_dcm(plan_file))
            logger.info(f"Loaded {len(plans)} plan files")
            return plans

    @cached_property
    def plan_summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan parameters, for all plans in the dataset.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan parameters, for all plans in the
            dataset.
        """
        return pl.concat(
            [
                plan.summary.with_columns(plan_index=pl.lit(key))
                for key, plan in enumerate(self.plans)
            ]
        ).select([pl.col("plan_index"), pl.all().exclude("plan_index")])

    @cached_property
    def plan_details(self) -> list[RadixactPlanDetails]:
        """Returns list containing treatment plan details.

        Returns
        -------
        list[RadixactPlanDetails]
            List of treatment plan details in the dataset.
        """
        plan_details = []
        plan_detail_files = self._get_filtered_series_list(
            "Plan Details", series="curr_path"
        )
        if len(plan_detail_files) == 0:
            logger.debug("Dataset contains no plan detail files")
            return plan_details
        else:
            for plan_detail_file in plan_detail_files:
                logger.debug(f"Loading {plan_detail_file}")
                plan_details.append(RadixactPlanDetails.from_xml(plan_detail_file))
            logger.info(f"Loaded {len(plan_details)} plan detail files")
            return plan_details

    @cached_property
    def plan_informations(self) -> list[RadixactPlanInformation]:
        """Returns list containing treatment plan information.

        Returns
        -------
        list[RadixactPlanInformation]
            List of treatment plan informations in the dataset.
        """
        plan_informations = []
        plan_information_files = self._get_filtered_series_list(
            "Plan Information", series="curr_path"
        )
        if len(plan_information_files) == 0:
            logger.debug("Dataset contains no plan information files")
            return plan_informations
        else:
            for plan_information_file in plan_information_files:
                logger.debug(f"Loading {plan_information_file}")
                plan_informations.append(
                    RadixactPlanInformation.from_xml(plan_information_file)
                )
            logger.info(f"Loaded {len(plan_informations)} plan information files")
            return plan_informations

    @cached_property
    def plan_settings(self) -> list[RadixactPlanSettings]:
        """Returns list containing treatment plan settings.

        Returns
        -------
        list[RadixactPlanSettings]
            List of treatment plan settings in the dataset.
        """
        plan_settings = []
        plan_setting_files = self._get_filtered_series_list(
            "Plan Settings", series="curr_path"
        )
        if len(plan_setting_files) == 0:
            logger.debug("Dataset contains no plan setting files")
            return plan_settings
        else:
            for plan_setting_file in plan_setting_files:
                logger.debug(f"Loading {plan_setting_file}")
                plan_settings.append(
                    RadixactPlanSettings.from_plan_settings(plan_setting_file)
                )
            logger.info(f"Loaded {len(plan_settings)} plan setting files")
            return plan_settings

    @cached_property
    def plan_settings_summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan settings, for all plan settings in dataset.

        Returns
        -------
        pl.DataFrame
            Summary of treatment plan settings, for all plan settings in dataset.
        """
        return pl.concat(
            [
                plan_settings.summary.with_columns(plan_settings_id=pl.lit(key))
                for key, plan_settings in enumerate(self.plan_settings)
            ]
        )

    @cached_property
    def plan_sinograms(self) -> list[RadixactSinogram]:
        """Returns list containing plan sinograms.

        Returns
        -------
        list[RadixactSinogram]
            List of treatment plan sinograms in the dataset.
        """
        # TODO Potentially include orig_path or curr_path.split("~")[-2] as fraction number
        # or include entire UID from orig_path or curr_path.
        plan_sinograms = []
        plan_sinogram_files = self._get_filtered_series_list(
            "Plan Sinogram", series="curr_path"
        )
        logger.debug(f"Dataset contains {len(plan_sinogram_files)} plan sinogram files")
        if len(plan_sinogram_files) == 0:
            return plan_sinograms
        else:
            for plan_sinogram_file in plan_sinogram_files:
                logger.debug(f"Loading {plan_sinogram_file}")
                plan_sinograms.append(RadixactSinogram.from_dplan(plan_sinogram_file))
            logger.info(f"Loaded {len(plan_sinograms)} plan sinogram files")
            return plan_sinograms

    @cached_property
    def radiation_records(self) -> list[pydicom.Dataset]:
        """Returns list containing fractional delivery radiation records.

        Returns
        -------
        list[pydicom.Dataset]
            List of radiation records in the dataset.
        """
        # TODO Potentially include ds.ReferenceRTPlanSequence[0].ReferenceSOPInstanceUID,
        # which matches ReferenceSOPInstanceUID in the RTPLAN, also ds.InstanceNumber,
        # which includes fraction_number as instance[:-2] and session as instance[-2:]
        radiation_records = []
        radiation_record_files = self._get_filtered_series_list(
            "Radiation Record", series="curr_path"
        )
        if len(radiation_record_files) == 0:
            logger.debug("Dataset contains no radiation record files")
            return radiation_records
        else:
            for radiation_record_file in radiation_record_files:
                logger.debug(f"Loading {radiation_record_file}")
                radiation_records.append(pydicom.dcmread(radiation_record_file))
            logger.info(f"Loaded {len(radiation_records)} radiation record files")
            return radiation_records

    @cached_property
    def structure_sets(self) -> list[pydicom.Dataset]:
        """Returns list containing structure sets.

        Returns
        -------
        list[pydicom.Dataset]
            List of structure sets in the dataset.
        """
        structure_sets = []
        structure_set_files = self._get_filtered_series_list(
            "Structure Set", series="curr_path"
        )
        if len(structure_set_files) == 0:
            logger.debug("Dataset contains no structure set files")
            return structure_sets
        else:
            for structure_set_file in structure_set_files:
                logger.debug(f"Loading {structure_set_file}")
                structure_sets.append(pydicom.dcmread(structure_set_file))
            logger.info(f"Loaded {len(structure_sets)} structure set files")
            return structure_sets

    @cached_property
    def telemetry_metrics(self) -> pl.DataFrame:
        results = []
        for session_index, telemetry_sinogram in enumerate(self.telemetry_sinograms):
            plan_sinogram = self._find_matching_planned_sinogram(telemetry_sinogram)
            if plan_sinogram is None:
                continue
            error_projections = plan_sinogram.detect_fluence_errors(telemetry_sinogram)
            fluence_variation_projections = plan_sinogram.detect_fluence_variations(
                telemetry_sinogram
            )
            error_transitions = np.diff(
                np.concatenate(([False], error_projections, [False])).astype(int)
            )
            error_lengths = (
                np.where(error_transitions == -1)[0]
                - np.where(error_transitions == 1)[0]
            )
            results.append(
                {
                    "session_index": int(session_index),
                    "num_fluence_variations": len(error_lengths),
                    "num_projections_with_fluence_variations": int(
                        np.sum(error_projections)
                    ),
                    "num_projections_with_adaptation": int(
                        np.sum(fluence_variation_projections)
                    ),
                    "longest_fluence_variation_duration": int(np.max(error_lengths))
                    if error_lengths.size > 0
                    else 0,
                    "longest_fluence_variation_duration_ms": int(np.max(error_lengths))
                    * telemetry_sinogram.projection_duration()
                    if error_lengths.size > 0
                    else 0.0,
                    "planned_zero_lot_projections": np.count_nonzero(
                        np.sum(plan_sinogram.fractional_leaf_open_times(), axis=1) == 0
                    ),
                    "telemetry_zero_lot_projections": np.count_nonzero(
                        np.sum(telemetry_sinogram.fractional_leaf_open_times(), axis=1)
                        == 0
                    ),
                }
            )
        return pl.DataFrame(results)

    @cached_property
    def telemetry_sinograms(self) -> list[RadixactSinogram]:
        """Returns list containing telemetry sinograms from fractional deliveries.

        Returns
        -------
        list[RadixactSinogram]
            List of telemetry sinograms in the dataset.
        """
        # TODO Potentially include orig_path or curr_path.split("-")[-3] as plan uid,
        # [-2] as fraction and [1] as fragment, or include entire UID from orig_path
        # or curr_path.
        telemetry_sinograms = []
        telemetry_sinogram_files = self._get_filtered_series_list(
            "Telemetry Sinogram", series="curr_path"
        )
        if len(telemetry_sinogram_files) == 0:
            logger.debug("Dataset contains no telemetry sinogram files")
            return telemetry_sinograms
        else:
            for telemetry_sinogram_file in telemetry_sinogram_files:
                telemetry_sinograms.append(
                    RadixactSinogram.from_dplan(telemetry_sinogram_file)
                )
            return telemetry_sinograms

    @cached_property
    def telemetry_timings(self) -> list[RadixactTiming]:
        """Returns list containing telemetry timings from fractional deliveries.

        Returns
        -------
        list[RadixactTiming]
            List of telemetry timings in the dataset.
        """
        telemetry_timings = []
        telemetry_timing_files = self._get_filtered_series_list(
            "Telemetry Timing", series="curr_path"
        )
        if len(telemetry_timing_files) == 0:
            return telemetry_timings
        else:
            for telemetry_timing_file in telemetry_timing_files:
                telemetry_timings.append(RadixactTiming.from_dat(telemetry_timing_file))
            return telemetry_timings

    # endregion

    # region Public methods

    def plot_motion_histogram(
        self,
        mode: str = "absolute",
        fig_size: tuple[float, float] = (12, 4),
        offset_lim: tuple[float, float] | None = None,
        offset_bin: float = 0.25,
        vector_lim: tuple[float, float] | None = None,
        vector_bin: float = 0.25,
        title: str | None = None,
    ) -> mpl.Figure:
        """Generate a histogram plot of distribution of offset values in each dimension.

        Parameters
        ----------
        mode : str, optional
            Indicates whether to use absolute or relative offset. Default is "absolute".
        fig_size : tuple[float, float], optional
            Size of the figure. Default is (12, 4) inches.
        offset_lim : tuple[float, float], optional
            Minimum and maximum bin boundaries for offset. Default is None, in which
            case the limits will be calculated with floor(min) and ceiling(max),
            respectively.
        offset_bin : float, optional
            Bin width. Default is 0.25 mm.
        vector_lim : tuple[float, float], optional
            Minimum and maximum bin boundaries for vector. Default is None, in which
            case the limits will be 0 and ceiling(max).
        vector_bin : float, optional
            Bin width, Default is 0.25 mm.
        title : str, optional
            Title for figure. Default is None.

        Returns
        -------
        mpl.Figure
            Histogram of target offset values in each dimension.
        """
        return self.motion.plot_motion_histogram(
            mode, fig_size, offset_lim, offset_bin, vector_lim, vector_bin, title
        )

    def plot_session_motions(
        self,
        parameters: list[str] | None = None,
        col_wrap: int = 4,
    ) -> mpl.Figure:
        """Prepares relational line plot for parameter distances, with each session as
        a unique facet.

        Parameters
        ----------
        parameters : list[str], optional
            The Synchrony paramters to include in the figure. Possible options include:
            "target_offset_x", "target_offset_y", "target_offset_z"
                Target offset in IEC-X, IEC-Y and/or IEC-Z dimensions.
            "target_offset_vector"
                Target offset in 3D space.
            "potential_diff"
                Potential difference or uncertainty in target position calculcation.
            "measured_diff"
                Measured delta or difference between observed and calculated offset.
            "rigid_body"
                Rigid body value describing deformation of fiducial distribution.
            Default is None, in which case target_offset_x, target_offset_y,
            target_offset_z, potential_diff, measured_diff and rigid_body are used.
        col_wrap : int, optional
            Number of columns in figure. Default is 4.

        Returns
        -------
        mpl.Figure
            Relational line plot for parameter distances.
        """
        if parameters is None:
            parameters = (
                [
                    "target_offset_x",
                    "target_offset_y",
                    "target_offset_z",
                    "target_offset_vector",
                    "potential_diff",
                    "measured_diff",
                    "rigid_body",
                ],
            )
        mapping = {
            "target_offset_x": "IEC-X target offset",
            "target_offset_y": "IEC-Y target offset",
            "target_offset_z": "IEC-Z target offset",
            "target_offset_vector": "Vector target offset",
            "potential_diff": "Potential difference",
            "measured_diff": "Measured delta",
            "rigid_body": "Rigid body",
        }
        select = [
            pl.col(col)
            for col in [
                "session_index",
                "delta_time",
            ]
            + parameters
            if col in self.motion._df
        ]
        unpivot_df = (
            self.motion._df.select(select)
            .unpivot(index=["session_index", "delta_time"])
            .select(
                [
                    (pl.col("session_index") + 1).alias("Session"),
                    pl.col("variable").replace(mapping).alias("Parameter"),
                    pl.col("value").alias("Distance (mm)"),
                    pl.col("delta_time").alias("Time (s)"),
                ]
            )
        )
        ax = sns.relplot(
            unpivot_df,
            y="Distance (mm)",
            x="Time (s)",
            col="Session",
            hue="Parameter",
            col_wrap=col_wrap,
            kind="line",
        )
        sns.move_legend(ax, "upper center", ncol=10, bbox_to_anchor=(0.5, 1.02))
        return ax.figure

    def plot_session_motions_fraction_less_than_threshold(
        self,
        offset_type: str = "target_offset_vector",
        threshold_step: float = 1,
        figsize=(12, 4),
    ) -> mpl.Figure:
        """Plots the fraction of session motion statistics less than thresholds,
        intended to provide insight on margin selection.

        Parameters
        ----------
        offset_type : str, optional
            The motion metric percentiles to evaluate against thresholds. Default is
            "target_offset_vector".
        threshold_step : float, optional
            Step width to use to define thresholds against which to evaluate offset
            percentiles, in mm. Default is 1.
        figsize : tuple, optional
            Figure size. Default is (12, 4), i.e., 12 by 4 inches.

        Returns
        -------
        mpl.Figure
            Figure showing fraction of session motion statistics less than thresholds.

        Notes
        -----
        This calculation is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        return self.motion.plot_session_fraction_less_than_threshold(
            offset_type, threshold_step, figsize
        )

    def save_compressed(
        self, path: str | os.PathLike, types: list[str] | None, overwrite: bool = True
    ) -> None:
        """Save compressed versions of dataset files in a directory.

        Parameters
        ----------
        path : str | os.PathLike
            Path to which the compressed files should be saved.
        types : list[str] | None, optional
            List of file types to save copies of. Default is None, in which case
            Detector Sinogram, Motion Data, Plan, Plan Details, Plan Information,
            Plan Sinogram, Record, Telemetry Sinogram and Telemetry Timing files
            will be copied.
        overwrite : bool, optional
            Flag to silently overwrite existing files. Default is True. If False,
            file copy will be silently skipped.
        """
        if types is None:
            types = [
                "Detector Sinogram",
                "Motion Data",
                # "Plan",
                # "Plan Details",
                # "Plan Information",
                # "Plan Settings",
                "Plan Sinogram",
                # "Record",
                "Telemetry Sinogram",
                # "Telemetry Timing",
            ]
        properties = {
            "Detector Sinogram": self.detector_sinograms,
            "Motion Data": self.motions,
            # "Plan": self.plans,
            # "Plan Details": self.plan_details,
            # "Plan Information": self.plan_informations,
            # "Plan Settings": self.plan_settings,
            "Plan Sinogram": self.plan_sinograms,
            # "Record": self.radiation_records,
            "Telemetry Sinogram": self.telemetry_sinograms,
            # "Telemetry Timing": self.telemetry_timings,
        }
        for type in types:
            for index, data in enumerate(properties[type]):
                dst_path = os.path.join(
                    path,
                    f"{self._DEFAULT_PREFIXES[type]}_{index + 1:03d}"
                    + f".{self._COMPRESSED_EXTENSIONS_EXTENSIONS[type]}",
                )
                if os.path.exists(dst_path):
                    if overwrite:
                        data.to_compressed(dst_path)
                else:
                    data.to_compressed(dst_path)

    def save_copy(
        self,
        path: str | os.PathLike,
        types: list[str] | None = None,
        overwrite: bool = True,
    ) -> None:
        """Save copies of dataset files in a directory.

        Parameters
        ----------
        path : str | os.PathLike
            Path to which the copied files should be saved.
        types : list[str], optional
            List of file types to save copies of. Default is None, in which case
            Detector Sinogram, Motion Data, Plan, Plan Details, Plan Information,
            Plan Sinogram, Record, Telemetry Sinogram and Telemetry Timing files
            will be copied.
        overwrite : bool, optional
            Flag to silently overwrite existing files. Default is True. If False,
            file copy will be silently skipped.
        """
        if types is None:
            types = [
                "Detector Sinogram",
                "Motion Data",
                "Plan",
                "Plan Details",
                "Plan Information",
                "Plan Settings",
                "Plan Sinogram",
                "Record",
                "Telemetry Sinogram",
                "Telemetry Timing",
            ]
        for type in types:
            for index, src_path in enumerate(self._get_filtered_series_list(type)):
                dst_path = os.path.join(
                    path,
                    f"{self._DEFAULT_PREFIXES[type]}_{index + 1:03d}"
                    + f".{self._DEFAULT_EXTENSIONS[type]}",
                )
                if os.path.exists(dst_path):
                    if overwrite:
                        shutil.copy(src_path, dst_path)
                else:
                    shutil.copy(src_path, dst_path)

    # endregion

    # region Private methods

    def _find_matching_planned_sinogram(self, telemetry_sinogram) -> RadixactSinogram:
        telemetry_sinogram_total_lots = np.sum(
            telemetry_sinogram.fractional_leaf_open_times(), axis=1
        )
        # TODO identify correct plan more effeciently
        # TODO handle incomplete delivery sessions
        matching_sinogram = None
        matching_sinogram_total_lot_difference = np.inf
        for plan_sinogram in self.plan_sinograms:
            if len(plan_sinogram) == len(telemetry_sinogram):
                plan_sinogram_total_lots = np.sum(
                    plan_sinogram.fractional_leaf_open_times(), axis=1
                )
                plan_sinogram_total_lot_difference = np.sum(
                    np.abs(plan_sinogram_total_lots - telemetry_sinogram_total_lots)
                )
                if (
                    plan_sinogram_total_lot_difference
                    < matching_sinogram_total_lot_difference
                ):
                    matching_sinogram = plan_sinogram
                    matching_sinogram_total_lot_difference = (
                        plan_sinogram_total_lot_difference
                    )
        return matching_sinogram

    def _get_filtered_series_list(self, type: str, series="curr_path") -> list[str]:
        """Returns a list of paths for files with specified type filter.

        Parameters
        ----------
        type : str
            The type (of file) to use as a filter.
        series : str, optional
            The series name, or column name, of the series to be returned as a list.
            Default is "curr_path".

        Returns
        -------
        list[str]
            List of file paths for filtered type.
        """
        # TODO verify series / column name is valid
        result = self._files.filter(pl.col("type") == type).get_column(series).to_list()
        return result

    # endregion

    # region Static methods

    @staticmethod
    def _extractor_data_files_dataframe(path: str, filter: bool = True) -> pl.DataFrame:
        """Read file system BO for all patient files contained in a directory exported
        using the Accuray Patient Data Extractor software.

        Parameters
        ----------
        path : str
            Path of directory contain patient data extractor files.
        filter : bool, optional
            Flag to remove unknown or redundant files from the resultant dataframe.
            Default is true.

        Returns
        -------
        pl.DataFrame
            DataFrame containing all dataset files contained in a specified path.
        """
        # instantiate result and schema
        result = []
        schema = [
            "key",
            "orig_path",
            "token",
            "curr_path",
            "status",
            "format",
            "length",
            "hash",
            "category",
            "insert_date",
            "insert_user_id",
            "last_insert_date",
            "storage_location_key",
        ]
        # find export data file manifest
        if os.path.exists(os.path.join(path, "PatientExportDataBO.xml")):
            tree = et.parse(os.path.join(path, "PatientExportDataBO.xml"))
        else:
            tree = et.parse(
                glob.glob(os.path.join(path, "PatientExportDataBO*.xml"))[0]
            )
        # iterate over files listed in manifest
        for file_child in tree.getroot().findall(
            "PatientFiles/FileSystemBOList/FileSystemBO"
        ):
            result.append(
                [
                    file_child.find("FileKey").text,
                    file_child.find("OrigFileName").text,
                    file_child.find("FileToken").text,
                    os.path.join(path, file_child.find("FileToken").text.split("$")[1]),
                    file_child.find("FileStatus").text,
                    file_child.find("FileFormatCde").text,
                    file_child.find("FileLength").text,
                    file_child.find("FileHash").text,
                    file_child.find("FileCategory").text,
                    file_child.find("InsertDte").text,
                    file_child.find("InsertUserIdKey").text,
                    file_child.find("LastInsertDte").text,
                    file_child.find("StorageLocationKey").text,
                ]
            )
        # prepare dataframe
        result = pl.DataFrame(data=result, schema=schema, orient="row")
        result = RadixactDataset._categorise_data_files_dataframe(
            result, series_name="orig_path"
        )
        # return filtered or unfiltered dataframe
        if filter:
            return result.filter(pl.col("type").str.contains("Unknown").not_())
        else:
            return result

    @staticmethod
    def _categorise_data_files_dataframe(
        df=pl.DataFrame, series_name: str = "curr_path"
    ) -> pl.DataFrame:
        """Categorises files in manifest dataframe with a "type" string, e.g.,
        "Plan", "Dose Distribution", "Structure Set", "Telemetry sinogram",
        to allow for simple filtering.

        Parameters
        ----------
        df : DataFrame
            Dataframe containing manifest of files defining a dataset.
        series_name : str, optional
            Name of series or column to be used for categorisation. Default is
            "curr_path".

        Returns
        -------
        pl.DataFrame
            DataFrame containing file paths categorised according to file type.
        """
        conditions = {
            "Detector Sinogram": pl.col(series_name).str.contains(
                "(?i)detector(byproj|fluence)"
            ),
            "Dose Distribution": pl.col(series_name).str.contains("(?i)rtdose"),
            "Motion Data": pl.col(series_name).str.contains("(?i)motiondata"),
            "Optimizer Dose Distribution": pl.col(series_name).str.contains(
                "(?i)optimizationdoseinfo_rtdose"
            ),
            "PreciseART Dose Distribution": pl.col(series_name).str.contains(
                "(?i)MIM_SERIES.*rtdose"
            ),
            "Plan": pl.col(series_name).str.contains("(?i)rtplan"),
            "Plan Details": pl.col(series_name).str.contains("(?i)plandetails"),
            "Plan Information": pl.col(series_name).str.contains(
                "(?i)(generalplan|planinfo)"
            ),
            "Plan Settings": pl.col(series_name).str.contains("(?i)plan.*settings"),
            "Plan Sinogram": pl.col(series_name).str.contains(
                "(?i)(finaldeliveryplan|planfluence)"
            ),
            "Record": pl.col(series_name).str.contains("(?i)rtrecord"),
            "Structure Set": pl.col(series_name).str.contains("(?i)rtss"),
            "Telemetry Sinogram": pl.col(series_name).str.contains("(?i)telemfluence"),
            "Telemetry Timing": pl.col(series_name).str.contains("(?i)telemtiming"),
        }
        expr = pl.lit("Unknown")
        for type, condition in conditions.items():
            expr = pl.when(condition).then(pl.lit(type)).otherwise(expr)
        return df.with_columns(expr.alias("type"))

    # end region
