# -*- coding: utf-8 -*-
"""RadixactDataset module.

This module provides functionality for processing of data from Radixact treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import glob
import logging
import os
import xml.etree.ElementTree as et
from functools import cached_property

import polars as pl
import pydicom

from pyradlib.radixact.motion import RadixactSynchronyMotion
from pyradlib.radixact.plan import RadixactPlanSettings
from pyradlib.radixact.sinogram import RadixactSinogram
from pyradlib.radixact.timing import RadixactTiming

logger = logging.getLogger(__name__)


class RadixactDataset:
    # region Constructors

    def __init__(self, label: str = None):
        self.label = label
        self.files = None
        """
        self.plans = {}
        self.plan_reports = {}
        self.plan_settings = {}
        self.fraction_records = {}
        self.fraction_motions = {}
        self.fraction_telemetry_sinograms = {}
        self.fraction_detector_sinograms = {}
        """

    @classmethod
    def from_data_extractor(cls, path):
        """

        Parameters
        ----------
        path : str
            Path to the patient data folder.

        Returns
        -------
        RadixactDataset

        Notes
        -----
        The specified path should correspond to a Patient Data Extractor folder,
        containing a PatientExportDataBO.xml file.
        """
        dataset = cls()
        dataset.files = cls._extractor_data_files_dataframe(path)
        logger.debug(f"Dataset contains {len(dataset.files)} files")
        return dataset

    @classmethod
    def from_delivery_analysis(cls, path):
        """

        Parameters
        ----------
        path : str
            Path containing the Delivery Analysis cache files to be read.

        Returns
        -------
        RadixactDataset

        Notes
        -----
        The specified directory should correspond to a Delivery Analysis cache
        folder, e.g., containing DetectorByProj*.det, TelemFluence_*.dplan, and
        FinalDeliveryPlan*.dplan files.
        """
        return cls.from_path_list(sorted(os.listdir(path)))

    @classmethod
    def from_path_list(cls, paths: list[str]):
        """
        Parameters
        ----------
        paths : list[str]
            Paths of the files to be included in the DataSet.

        Returns
        -------
        RadixactDataset
        """
        dataset = cls()
        dataset.files = cls._categorise_data_files_dataframe(
            pl.DataFrame({"curr_path": paths}), series_name="curr_path"
        )
        return dataset

    # endregion

    # region Magic methods

    # endregion

    # region Properties

    @cached_property
    def detector_sinograms(self):
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
    def dose_distributions(self):
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
        if len(self.motions) == 0:
            return None
        else:
            result = sum(self.motions)
            logger.info(f"Combined {len(self.motions)} motion data files")
            return result

    @cached_property
    def motions(self) -> list[RadixactSynchronyMotion]:
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
                motions.append(RadixactSynchronyMotion.from_xml(motion_file))
            logger.info(f"Loaded {len(motions)} motion data files")
            return motions

    @cached_property
    def plans(self) -> list[pydicom.Dataset]:
        # TODO Potentially extract ds.BeamSequence[0].DeviceSerialNumber and
        # ds.RTPlanLabel to be used as identifiers in a dictionary.
        plans = []
        plan_files = self._get_filtered_series_list(
            "Treatment Plan", series="curr_path"
        )
        if len(plan_files) == 0:
            logger.debug("Dataset contains no plan files")
            return plans
        else:
            for plan_file in plan_files:
                logger.debug(f"Loading {plan_file}")
                plans.append(pydicom.dcmread(plan_file))
            logger.info(f"Loaded {len(plans)} plan files")
            return plans

    @cached_property
    def plan_details(self):
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
                plan_details.append(
                    RadixactPlanSettings.from_plan_settings(plan_detail_file)
                )
            logger.info(f"Loaded {len(plan_details)} plan detail files")
            return plan_details

    @cached_property
    def plan_settings(self):
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
    def plan_sinograms(self) -> RadixactSinogram:
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
    def radiation_records(self):
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
    def structure_sets(self):
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
    def telemetry_sinograms(self) -> list[RadixactSinogram]:
        # TODO Potentially include orig_path or curr_path.split("-")[-3] as plan uid,
        # [-2] as fraction and [1] as fragment, or include entire UID from orig_path
        # or curr_path.
        telemetry_sinograms = []
        telemetry_sinogram_files = self._get_filtered_series_list(
            "Telemetry Sinogram", series="curr_path"
        )
        telemetry_sinogram_orig_paths = self._get_filtered_series_list(
            "Telemetry Sinogram", series="orig_path"
        )
        if len(telemetry_sinogram_files) == 0:
            logger.debug("Dataset contains no telemetry sinogram files")
            return telemetry_sinograms
        else:
            if len(telemetry_sinogram_files) == len(self.telemetry_timings):
                telemetry_timing_dict = {
                    timing.uid: timing for timing in self.telemetry_timings
                }
            for (
                telemetry_sinogram_file,
                telemetry_sinogram_orig_path,
            ) in zip(
                telemetry_sinogram_files,
                telemetry_sinogram_orig_paths,
            ):
                uid = (
                    telemetry_sinogram_orig_path.split("\\")[-1]
                    .replace("TelemFluence_", "")
                    .replace(".dplan", "")
                )
                if len(telemetry_sinogram_files) == len(self.telemetry_timings):
                    logger.debug(f"Loading {telemetry_sinogram_file}, applying timing")
                    telemetry_sinograms.append(
                        RadixactSinogram.from_dplan(
                            telemetry_sinogram_file, uid
                        ).add_timestamps(telemetry_timing_dict[uid])
                    )
                else:
                    logger.debug(f"Loading {telemetry_sinogram_file}")
                    telemetry_sinograms.append(
                        RadixactSinogram.from_dplan(telemetry_sinogram_file, uid)
                    )
            return telemetry_sinograms

    @cached_property
    def telemetry_timings(self) -> list[RadixactTiming]:
        telemetry_timings = []
        telemetry_timing_files = self._get_filtered_series_list(
            "Telemetry Timing", series="curr_path"
        )
        telemetry_timing_orig_paths = self._get_filtered_series_list(
            "Telemetry Timing", series="orig_path"
        )
        if len(telemetry_timing_files) == 0:
            return telemetry_timings
        else:
            for telemetry_timing_file, telemetry_timing_orig_path in zip(
                telemetry_timing_files, telemetry_timing_orig_paths
            ):
                uid = (
                    telemetry_timing_orig_path.split("\\")[-1]
                    .replace("TelemTiming_", "")
                    .replace(".dat", "")
                )
                logger.debug(f"Loading {telemetry_timing_file}")
                telemetry_timings.append(
                    RadixactTiming.from_dat(telemetry_timing_file, uid)
                )
            logger.info(f"Loaded {len(telemetry_timings)} telemetry timing files")
            return telemetry_timings

    # endregion

    # region Public methods

    # endregion

    # region Private methods

    def _get_filtered_series_list(self, type: str, series="curr_path") -> list[str]:
        """
        Parameters
        ----------
        type : str
            The type (of file) to use as a filter.
        series : str, optional
            The series name, or column name, of the series to be returned as a list.
        """
        # TODO verify series / column name is valid
        result = self.files.filter(pl.col("type") == type).get_column(series).to_list()
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
        "Treatment Plan", "Dose Distribution", "Structure Set", "Telemetry sinogram",
        to allow for simple filtering.

        Parameters
        ----------
        df : DataFrame
            Dataframe containing manifest of files defining a dataset.
        series_name : str, optional
            Name of series or column to be used for categorisation. Default is
            "curr_path".
        """
        return df.with_columns(
            pl.when(pl.col(series_name).str.contains("(?i)detectorbyproj"))
            .then(pl.lit("Detector Sinogram"))
            .when(pl.col(series_name).str.contains("(?i)finaldoseinfo_rtdose"))
            .then(pl.lit("Dose Distribution"))
            .when(pl.col(series_name).str.contains("(?i)motiondata"))
            .then(pl.lit("Motion Data"))
            .when(pl.col(series_name).str.contains("(?i)rtplan"))
            .then(pl.lit("Plan"))
            .when(pl.col(series_name).str.contains("(?i)tomoplandetails"))
            .then(pl.lit("Plan Details"))
            .when(pl.col(series_name).str.contains("(?i)generalplan"))
            .then(pl.lit("Plan Information"))
            .when(pl.col(series_name).str.contains("(?i)plan.settings"))
            .then(pl.lit("Plan Settings"))
            .when(pl.col(series_name).str.contains("(?i)finaldeliveryplan"))
            .then(pl.lit("Plan Sinogram"))
            .when(pl.col(series_name).str.contains("(?i)rtrecord"))
            .then(pl.lit("Radiation Record"))
            .when(pl.col(series_name).str.contains("(?i)rtss"))
            .then(pl.lit("Structure Set"))
            .when(pl.col(series_name).str.contains("(?i)telemfluence"))
            .then(pl.lit("Telemetry Sinogram"))
            .when(pl.col(series_name).str.contains("(?i)telemtiming"))
            .then(pl.lit("Telemetry Timing"))
            .otherwise(pl.lit("Unknown"))
            .alias("type")
        )

    # end region
